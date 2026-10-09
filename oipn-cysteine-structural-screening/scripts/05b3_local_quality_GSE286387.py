#!/usr/bin/env python3
"""Offline local confidence gate; assembly/chemical/final structural approval remain pending."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import importlib.util
import json
import math
from pathlib import Path
import platform
import shutil
import sys
import tempfile

SPEC = importlib.util.spec_from_file_location('mapping', Path(__file__).with_name('05b2_map_cysteines_GSE286387.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
I = M.I
VERSION = '1.0.1'
ROUNDING_TOLERANCE = 0.011  # CIF and confidence JSON may independently round to 2 decimals.
CYS_ATOMS = {'N', 'CA', 'C', 'O', 'CB', 'SG'}


class QualityHold(ValueError):
    pass


def need(condition, reason):
    if not condition:
        raise QualityHold(reason)


def numeric(value):
    if isinstance(value, bool):
        return None
    return M.number(value)


def load_input(project):
    path = project / 'config/GSE286387_step05b3_input.json'
    b = json.loads(path.read_text())
    I.require(b['status'] == 'STEP05B2_MAPPING_ACCEPTED_FOR_LOCAL_CONFIDENCE' and
              b['protocol_version'] == '0.1.0', 'Unaccepted local-quality input')
    inputs = [path]
    for name, expected in b['input_sha256'].items():
        f = I.safe_path(project, name)
        I.require(f.is_file() and M.file_sha(f) == expected, 'Frozen quality input changed: ' + name)
        inputs.append(f)
    raw_binding, seqs, candidates, links, ledger, archive, upstream = M.load_input(project)
    inputs += upstream
    snapshot = I.safe_path(project, b['mapping_snapshot_directory'])
    mapped = I.read_csv(snapshot / 'site_structure_mapping.csv')
    audit = I.read_csv(snapshot / 'site_step_audit.csv')
    summary = json.loads((snapshot / 'summary.json').read_text())
    I.require((snapshot / 'SUCCESS.txt').read_text().strip() == summary['status'] ==
              'STEP05B2_GENERATED_REVIEW_PENDING', 'Mapping generation not complete')
    I.require(len(mapped) == summary['mapping_records'] == 25026 and len(audit) == 10799,
              'Mapping accounting changed')
    for row in I.read_csv(snapshot / 'input_checksums.csv'):
        f = I.safe_path(project, row['input'])
        I.require(M.file_sha(f) == row['sha256'], 'Mapping provenance changed: ' + row['input'])
        inputs.append(f)
    expected = M.site_audit(seqs, mapped, {c['protein_accession'] for c in candidates if c['candidate_status'] == 'held'})
    I.require(audit == [{k: str(v) for k, v in r.items()} for r in expected], 'Site mapping audit inconsistent')
    raw = {r['file_id']: r for r in ledger}
    for row in mapped:
        I.require(row['structure_file_id'] in raw and row['structure_sha256'] == raw[row['structure_file_id']]['sha256']
                  and row['local_file'] == raw[row['structure_file_id']]['local_file'], 'Mapping/raw provenance mismatch')
        c = candidates[int(row['candidate_row']) - 1]
        I.require(c['protein_accession'] == row['protein_accession'] and c['structure_id'] == row['structure_id']
                  and c['source'] == row['source'] and c['candidate_status'] == 'candidate', 'Mapping candidate mismatch')
    parameters = project / 'config/analysis_parameters.json'
    fixed = json.loads(parameters.read_text())['structure']
    I.require(fixed['local_radius_angstrom'] == 6.0 and fixed['predicted_target_plddt_minimum'] == 90
              and fixed['predicted_neighbor_plddt_minimum'] == 70 and fixed['predicted_neighbor_plddt_strict'] == 90,
              'Frozen local-confidence thresholds changed')
    inputs += [parameters, project / 'requirements-structure.txt', Path(__file__).resolve()]
    return b, seqs, candidates, links, ledger, archive, mapped, audit, fixed, list(dict.fromkeys(inputs))


def confidence_values(data, length):
    need(isinstance(data, dict), 'CONFIDENCE_SCHEMA_UNRESOLVED')
    nums, scores = data.get('residueNumber'), data.get('confidenceScore')
    need(isinstance(nums, list) and nums == list(range(1, length + 1)) and
         all(type(n) is int for n in nums) and isinstance(scores, list) and len(scores) == length,
         'CONFIDENCE_RESIDUE_NUMBERING_UNRESOLVED')
    values = [numeric(x) for x in scores]
    need(all(type(x) in (int, float) for x in scores) and
         all(x is not None and 0 <= x <= 100 for x in values), 'CONFIDENCE_VALUES_UNRESOLVED')
    return dict(zip(nums, values))


def pae_matrix(data, length):
    need(isinstance(data, list) and len(data) == 1 and isinstance(data[0], dict), 'PAE_SCHEMA_UNRESOLVED')
    matrix = data[0].get('predicted_aligned_error')
    maximum = numeric(data[0].get('max_predicted_aligned_error'))
    need(isinstance(matrix, list) and len(matrix) == length and maximum is not None and maximum >= 0,
         'PAE_DIMENSIONS_UNRESOLVED')
    for row in matrix:
        need(isinstance(row, list) and len(row) == length, 'PAE_DIMENSIONS_UNRESOLVED')
        # AFDB's integer-valued matrix rounds entries (e.g. 32) while max PAE may be 31.75.
        need(all(numeric(v) is not None and 0 <= numeric(v) <= maximum +
                 (0.5 if type(v) is int else ROUNDING_TOLERANCE) for v in row),
             'PAE_VALUES_UNRESOLVED')
    return matrix


def atom_index(block, label, model, confidence):
    """Keep the complete modeled monomer; no altloc, atom or missing-position repair."""
    atoms, by_residue, cells = [], defaultdict(list), defaultdict(list)
    for a in M.records(block, '_atom_site.'):
        need(a.get('label_asym_id') == label and (a.get('pdbx_PDB_model_num') or '1') == model,
             'PREDICTED_MODEL_NOT_SINGLE_MONOMER')
        need(a.get('label_seq_id', '').isdigit(), 'PREDICTED_NONPOLYMER_CONTEXT_UNRESOLVED')
        n = int(a['label_seq_id'])
        need(n in confidence, 'PREDICTED_ATOM_NUMBERING_UNRESOLVED')
        need(not a.get('label_alt_id') and numeric(a.get('occupancy')) == 1,
             'PREDICTED_ALTERNATE_OR_OCCUPANCY_CONTEXT_UNRESOLVED')
        xyz = tuple(numeric(a.get(k)) for k in ('Cartn_x', 'Cartn_y', 'Cartn_z'))
        need(all(v is not None for v in xyz), 'PREDICTED_COORDINATES_UNRESOLVED')
        score = numeric(a.get('B_iso_or_equiv'))
        need(score is not None and abs(score - confidence[n]) <= ROUNDING_TOLERANCE,
             'CIF_JSON_CONFIDENCE_DISAGREEMENT')
        element = a.get('type_symbol', '').upper()
        need(bool(element), 'ATOM_ELEMENT_UNRESOLVED')
        atom = dict(a, seq=n, xyz=xyz)
        atoms.append(atom)
        by_residue[n].append(atom)
        if element not in ('H', 'D'):
            cells[tuple(math.floor(v / 6.0) for v in xyz)].append(atom)
    need(set(by_residue) == set(confidence), 'PREDICTED_POLYMER_RESIDUES_UNOBSERVED')
    for residue in by_residue.values():
        need(len({a['label_atom_id'] for a in residue}) == len(residue), 'PREDICTED_DUPLICATE_ATOM_NAMES')
        need(len({a['label_comp_id'] for a in residue}) == 1, 'PREDICTED_MICROHETEROGENEITY')
    return atoms, by_residue, cells


def neighbor_atoms(cells, point, target, radius):
    center = tuple(math.floor(v / 6.0) for v in point)
    found = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                for a in cells.get((center[0] + dx, center[1] + dy, center[2] + dz), []):
                    if a['seq'] != target and sum((x - y) ** 2 for x, y in zip(a['xyz'], point)) <= radius ** 2:
                        found.append(a)
    return found


def initial_row(mapping):
    keep = ('site_id', 'protein_accession', 'canonical_cys_position', 'candidate_row', 'source', 'structure_id',
            'label_asym_id', 'model_id', 'label_seq_id', 'structure_file_id', 'structure_sha256')
    return dict({k: mapping[k] for k in keep}, mapping_status=mapping['mapping_status'],
                local_gate_status='held', reason_code='', target_plddt='', minimum_neighbor_plddt='',
                neighbor_residue_count='', neighbor_atom_count='', neighbor_canonical_positions='',
                neighbor_90_sensitivity='', coherent_target_atoms='', confidence_file_id='', confidence_sha256='',
                pae_file_id='', pae_sha256='', pae_status='not_assessed', max_target_neighbor_pae_bidirectional='',
                full_canonical_model='', fragment_boundary_review='',
                neighbor_atom_completeness='pending_context_review',
                assembly_context='pending_review', mutation_context='pending_review',
                structural_eligibility='held', protocol_version='0.1.0')


def evaluate_site(mapping, candidate, seq, confidence, matrix, by_residue, cells, fixed):
    row = initial_row(mapping)
    n = int(mapping['label_seq_id'])
    row['target_plddt'] = confidence[n]
    target = by_residue[n]
    need(all(a['label_comp_id'] == 'CYS' for a in target), 'TARGET_COMPONENT_UNRESOLVED')
    names = {a['label_atom_id'] for a in target if a['type_symbol'].upper() not in ('H', 'D')}
    row['coherent_target_atoms'] = ';'.join(sorted(names))
    if not CYS_ATOMS <= names:
        row['reason_code'] = 'COHERENT_TARGET_HEAVY_ATOMS_MISSING'
        return row
    sulfur = [a for a in target if a['label_atom_id'] == 'SG' and a['type_symbol'].upper() == 'S']
    need(len(sulfur) == 1, 'TARGET_SG_ATOM_UNRESOLVED')
    sg = sulfur[0]
    I.require(sg['id'] == mapping['selected_SG_atom_id'] and all(abs(x - float(mapping[k])) < 1e-6
         for x, k in zip(sg['xyz'], ('SG_x', 'SG_y', 'SG_z'))), 'MAPPED_SG_COORDINATES_DISAGREE')
    neighbors = neighbor_atoms(cells, sg['xyz'], n, fixed['local_radius_angstrom'])
    residues = sorted({a['seq'] for a in neighbors})
    start = int(candidate['sequence_start'])
    end = int(candidate['sequence_end'])
    I.require(start + n - 1 == int(mapping['canonical_cys_position']), 'MAPPED_CANONICAL_POSITION_DISAGREES')
    row.update(target_plddt=confidence[n], minimum_neighbor_plddt=min((confidence[i] for i in residues), default=''),
               neighbor_residue_count=len(residues), neighbor_atom_count=len(neighbors),
               neighbor_canonical_positions=';'.join(str(start + i - 1) for i in residues),
               neighbor_90_sensitivity=all(confidence[i] >= fixed['predicted_neighbor_plddt_strict'] for i in residues),
               full_canonical_model=start == 1 and end == len(seq), fragment_boundary_review=not (start == 1 and end == len(seq)))
    if matrix is not None:
        row.update(pae_status='parsed_descriptive_only', max_target_neighbor_pae_bidirectional=
                   max((max(matrix[n - 1][i - 1], matrix[i - 1][n - 1]) for i in residues), default=''))
    if confidence[n] < fixed['predicted_target_plddt_minimum']:
        row.update(local_gate_status='excluded', reason_code='TARGET_PLDDT_BELOW_90')
    elif any(confidence[i] < fixed['predicted_neighbor_plddt_minimum'] for i in residues):
        row.update(local_gate_status='excluded', reason_code='NEIGHBOR_PLDDT_BELOW_70')
    else:
        row.update(local_gate_status='pass', reason_code='LOCAL_CONFIDENCE_PASS_CONTEXT_PENDING')
    return row


def evaluate_model(block, candidate, seq, mapping, confidence_data, pae_data, fixed):
    """Expected scientific ambiguity holds the model; parser/file corruption fails the run."""
    try:
        entities, bad, types, asym = M.polymer_index(block)
        labels = {r['label_asym_id'] for r in mapping}
        models = {r['model_id'] for r in mapping}
        need(len(labels) == len(models) == 1 and '' not in labels and '' not in models,
             'PREDICTED_CHAIN_MODEL_CONTEXT_UNRESOLVED')
        label, model = next(iter(labels)), next(iter(models))
        entity = asym[label]
        need(entity not in bad and types.get(entity) in ('polypeptide(L)', 'polypeptide(D)'),
             'PREDICTED_POLYMER_CONTEXT_UNRESOLVED')
        positions = entities[entity]
        length = len(positions)
        need(sorted(positions) == list(range(1, length + 1)), 'PREDICTED_POLYMER_NUMBERING_UNRESOLVED')
        start, end = int(candidate['sequence_start']), int(candidate['sequence_end'])
        need(''.join(M.aa(positions[i]) for i in range(1, length + 1)) == seq[start - 1:end]
             and length == end - start + 1, 'PREDICTED_SEQUENCE_CONTEXT_DISAGREEMENT')
        scores = confidence_values(confidence_data, length)
        _, by_residue, cells = atom_index(block, label, model, scores)
        pae_reason = ''
        try:
            matrix = pae_matrix(pae_data, length)
        except QualityHold as e:
            matrix = None
            pae_reason = str(e)
        out = []
        for original in mapping:
            try:
                row = evaluate_site(original, candidate, seq, scores, matrix, by_residue, cells, fixed)
            except QualityHold as e:
                row = initial_row(original)
                row['reason_code'] = str(e)
            if pae_reason:
                row['pae_status'] = pae_reason  # Never replaced by zero; not a new confidence cutoff.
            out.append(row)
        return out
    except (QualityHold, M.MappingHold) as e:
        out = [initial_row(r) for r in mapping]
        for r in out:
            r['reason_code'] = str(e)
        return out


def upstream_reviews(sequences, candidates, mapped):
    """Keep unresolved offered metadata and per-option mapping evidence visible downstream."""
    mapping_reviews = Counter(r['site_id'] for r in mapped if r['mapping_status'] == 'held')
    accessions = {c['protein_accession'] for c in candidates if c['candidate_status'] == 'held'}
    metadata_reviews = {f'{accession}:C{n}' for accession in accessions
                        for n, residue in enumerate(sequences[accession], 1) if residue == 'C'}
    return mapping_reviews, metadata_reviews


def aggregate_sites(audit, rows, mapping_reviews=None, metadata_reviews=None):
    mapping_reviews = mapping_reviews or {}
    metadata_reviews = metadata_reviews or set()
    groups = defaultdict(list)
    for r in rows:
        groups[r['site_id']].append(r)
    out = []
    for site in audit:
        options = groups[site['object_id']]
        passing = sum(r['local_gate_status'] == 'pass' for r in options)
        if passing:
            status, reason = 'held', 'LOCAL_CONFIDENCE_OPTIONS_CONTEXT_PENDING'
        elif (any(r['local_gate_status'] == 'held' for r in options) or
              site['reason_code'] == 'MAPPING_EVIDENCE_UNRESOLVED' or
              mapping_reviews.get(site['object_id'], 0) or site['object_id'] in metadata_reviews):
            status, reason = 'held', 'LOCAL_OR_MAPPING_REVIEW_PENDING'
        elif options:
            status, reason = 'excluded', 'ALL_MAPPED_OPTIONS_FAIL_LOCAL_CONFIDENCE'
        else:
            status, reason = 'unassessable', site['reason_code']
        out.append(dict(step_id='05b3', object_level='site', object_id=site['object_id'],
                        protein_accession=site['protein_accession'], canonical_cys_position=site['canonical_cys_position'],
                        status=status, reason_code=reason, assessed_mapping_options=len(options),
                        local_confidence_pass_options=passing,
                        unresolved_mapping_options=mapping_reviews.get(site['object_id'], 0),
                        unresolved_metadata_review=site['object_id'] in metadata_reviews,
                        structural_eligibility='held' if status == 'held' else status,
                        protocol_version='0.1.0'))
    return out


def run(project):
    out = project / 'results/05b3_local_quality'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'SUCCESS.txt').unlink(missing_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='quality_stage_', dir=out))
    try:
        g = M.gemmi_module()
        b, seqs, candidates, links, ledger, archive, mapped, audit, fixed, inputs = load_input(project)
        print('[1/3] Rechecking all local raw hashes and accepted mapping provenance...', flush=True)
        verified = M.verify_raw(project, ledger, archive)
        raw = {r['file_id']: r for r in ledger}
        groups, rows = defaultdict(list), []
        for r in mapped:
            if r['mapping_status'] != 'mapped':
                continue
            if r['source'] == 'AlphaFold_DB':
                groups[int(r['candidate_row'])].append(r)
            else:
                q = initial_row(r)
                q['reason_code'] = 'EXPERIMENTAL_LOCAL_VALIDATION_AND_CONTEXT_PENDING'
                rows.append(q)
        print('[2/3] Checking predicted confidence JSON/CIF and SG neighborhoods...', flush=True)
        model_audit = []
        for number, (index, options) in enumerate(sorted(groups.items()), 1):
            candidate = candidates[index - 1]
            files = [raw[fid] for fid in links[index - 1]['file_ids'].split(';')]
            by_kind = {r['kind']: r for r in files}
            I.require(len(by_kind) == len(files) == 3 and set(by_kind) ==
                      {'alphafold_mmcif', 'alphafold_confidence', 'alphafold_pae'}, 'Predicted file association unresolved')
            coordinate = by_kind['alphafold_mmcif']
            block = g.cif.read_file(str(I.safe_path(project, coordinate['local_file']))).sole_block()
            I.require(block.name == candidate['structure_id'], 'Predicted coordinate identity changed')
            conf, pae = by_kind['alphafold_confidence'], by_kind['alphafold_pae']
            evaluated = evaluate_model(block, candidate, seqs[candidate['protein_accession']], options,
                                       json.loads(I.safe_path(project, conf['local_file']).read_text()),
                                       json.loads(I.safe_path(project, pae['local_file']).read_text()), fixed)
            for r in evaluated:
                r.update(confidence_file_id=conf['file_id'], confidence_sha256=conf['sha256'],
                         pae_file_id=pae['file_id'], pae_sha256=pae['sha256'])
            rows += evaluated
            model_audit.append(dict(candidate_row=index, protein_accession=candidate['protein_accession'],
                                    structure_id=candidate['structure_id'], evaluated_sites=len(evaluated),
                                    local_gate_status_counts=json.dumps(dict(Counter(r['local_gate_status'] for r in evaluated)), sort_keys=True),
                                    reason_counts=json.dumps(dict(Counter(r['reason_code'] for r in evaluated)), sort_keys=True),
                                    structural_eligibility='held', protocol_version='0.1.0'))
            del block
            if number % 25 == 0 or number == len(groups):
                print(f'  Predicted models checked: {number}/{len(groups)}', flush=True)
        rows.sort(key=lambda r: (r['site_id'], int(r['candidate_row']), r['label_asym_id'], r['model_id']))
        I.require(len(rows) == sum(r['mapping_status'] == 'mapped' for r in mapped), 'Mapped-option accounting lost')
        mapping_reviews, metadata_reviews = upstream_reviews(seqs, candidates, mapped)
        sites = aggregate_sites(audit, rows, mapping_reviews, metadata_reviews)
        proteins = []
        for accession in sorted(seqs):
            subset = [r for r in sites if r['protein_accession'] == accession]
            proteins.append(dict(protein_accession=accession, input_cys_sites=len(subset),
                                 sites_with_local_confidence_options=sum(r['local_confidence_pass_options'] > 0 for r in subset),
                                 structural_eligibility='held', protocol_version='0.1.0'))
        for name, data in [('site_local_quality.csv', rows), ('predicted_model_audit.csv', model_audit),
                           ('site_step_audit.csv', sites), ('protein_quality_coverage.csv', proteins),
                           ('raw_file_checksums.csv', verified)]:
            I.require(bool(data), 'Unexpected empty quality output')
            I.write_csv(stage / name, data, list(data[0]))
        I.write_csv(stage / 'input_checksums.csv', [dict(input=str(f.relative_to(project)), sha256=M.file_sha(f)) for f in inputs],
                    ['input', 'sha256'])
        summary = dict(status='STEP05B3_GENERATED_REVIEW_PENDING', implementation_version=VERSION, protocol_version='0.1.0',
                       input_proteins=len(proteins), input_sites=len(sites), verified_raw_files=len(verified),
                       predicted_models_checked=len(groups), mapped_options=len(rows),
                       local_gate_status_counts=dict(Counter(r['local_gate_status'] for r in rows)),
                       reason_counts=dict(Counter(r['reason_code'] for r in rows)),
                       site_status_counts=dict(Counter(r['status'] for r in sites)),
                       sites_with_local_confidence_options=sum(r['local_confidence_pass_options'] > 0 for r in sites),
                       proteins_with_local_confidence_options=sum(r['sites_with_local_confidence_options'] > 0 for r in proteins),
                       final_structural_passes=0, python=platform.python_version(), gemmi=g.__version__,
                       generated_utc=datetime.now(timezone.utc).isoformat())
        (stage / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        report = '\n'.join(['# Step 05b3 — local predicted confidence and target completeness', '',
                            'Status: ' + summary['status'], '',
                            f'- Input proteins: {len(proteins)}; canonical Cys sites: {len(sites)}',
                            f'- Predicted models checked: {len(groups)}; mapped options retained: {len(rows)}',
                            '- Local gate statuses: ' + str(summary['local_gate_status_counts']),
                            f'- Sites with local-confidence options: {summary["sites_with_local_confidence_options"]}',
                            f'- Proteins with local-confidence options: {summary["proteins_with_local_confidence_options"]}', '',
                            '- Canonical site audit statuses: ' + str(summary['site_status_counts']), '',
                            'Frozen gate: target pLDDT >=90; every modeled non-H/D neighbor residue within 6 A of SG >=70.',
                            'The target residue is excluded from neighbor counts. JSON/CIF agreement and coherent CYS heavy atoms are required.',
                            'Neighbor >=90 is a sensitivity flag only. PAE is descriptive; no new PAE threshold is introduced.',
                            'Experimental mapped options remain held for manual local validation and assembly/mutation context.',
                            'Every original canonical site is retained, including mapping-unresolved and no-structure cases.',
                            'Unresolved offered metadata/mapping options remain held even if the mapped predicted option fails.',
                            'Local gate pass is not final structural approval: assembly, native/fragment context and chemical state remain pending.',
                            'No model selection, SASA, pKa or oxidation score is generated.', ''])
        (stage / 'quality_report.md').write_text(report)
        for f in stage.iterdir():
            f.replace(out / f.name)
        (out / 'FAILURE.txt').unlink(missing_ok=True)
        (out / 'SUCCESS.txt').write_text(summary['status'] + '\n')
        print('[3/3] Outputs: ' + str(out), flush=True)
        print(report)
        return summary
    except Exception as e:
        (out / 'FAILURE.txt').write_text(str(e) + '\n')
        raise
    finally:
        shutil.rmtree(stage, ignore_errors=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        run(args.project.resolve())
    except Exception as e:
        print('STEP05B3 FAILED: ' + str(e), file=sys.stderr)
        sys.exit(1)
