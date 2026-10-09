#!/usr/bin/env python3
"""Collect structural/experimental context evidence; do not select or approve a structure."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import importlib.util
import json
import math
from pathlib import Path
import platform
import shutil
import sys
import tempfile


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


Q = module('quality', '05b3_local_quality_GSE286387.py')
V = module('validation', '05b4_validation_support.py')
M, I = Q.M, Q.I
VERSION = '1.0.0'
BACKBONE = {'N', 'CA', 'C', 'O'}
SIDECHAINS = {
    'GLY': '', 'ALA': 'CB', 'VAL': 'CB CG1 CG2', 'LEU': 'CB CG CD1 CD2', 'ILE': 'CB CG1 CG2 CD1',
    'SER': 'CB OG', 'THR': 'CB OG1 CG2', 'CYS': 'CB SG', 'MET': 'CB CG SD CE',
    'PRO': 'CB CG CD', 'ASP': 'CB CG OD1 OD2', 'ASN': 'CB CG OD1 ND2',
    'GLU': 'CB CG CD OE1 OE2', 'GLN': 'CB CG CD OE1 NE2', 'LYS': 'CB CG CD CE NZ',
    'ARG': 'CB CG CD NE CZ NH1 NH2', 'HIS': 'CB CG ND1 CD2 CE1 NE2',
    'PHE': 'CB CG CD1 CD2 CE1 CE2 CZ', 'TYR': 'CB CG CD1 CD2 CE1 CE2 CZ OH',
    'TRP': 'CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2', 'SEC': 'CB SE',
}
ANNOTATIONS = ('_entity.', '_entity_src_gen.', '_entity_src_nat.', '_pdbx_entity_src_syn.',
               '_pdbx_struct_assembly.', '_pdbx_struct_assembly_gen.', '_pdbx_struct_oper_list.',
               '_struct_ref_seq_dif.', '_pdbx_unobs_or_zero_occ_residues.', '_pdbx_unobs_or_zero_occ_atoms.',
               '_struct_conn.', '_pdbx_audit_revision_history.')


def load_input(project):
    path = project / 'config/GSE286387_step05b4_input.json'
    b = json.loads(path.read_text())
    I.require(b['status'] == 'STEP05B3_LOCAL_CONFIDENCE_ACCEPTED_FOR_CONTEXT_REVIEW' and
              b['protocol_version'] == '0.1.0', 'Unaccepted context input')
    inputs = [path]
    for name, expected in b['input_sha256'].items():
        f = I.safe_path(project, name)
        I.require(f.is_file() and M.file_sha(f) == expected, 'Frozen context input changed: ' + name)
        inputs.append(f)
    _, sequences, candidates, links, ledger, archive, mapped, prior, fixed, upstream = Q.load_input(project)
    snapshot = I.safe_path(project, b['local_quality_snapshot_directory'])
    summary = json.loads((snapshot / 'summary.json').read_text())
    I.require((snapshot / 'SUCCESS.txt').read_text().strip() == summary['status'] ==
              'STEP05B3_GENERATED_REVIEW_PENDING' and summary['implementation_version'] == '1.0.1',
              'Corrected local-quality generation missing')
    for r in I.read_csv(snapshot / 'input_checksums.csv'):
        f = I.safe_path(project, r['input'])
        I.require(M.file_sha(f) == r['sha256'], 'Local-quality provenance changed: ' + r['input'])
        inputs.append(f)
    quality = I.read_csv(snapshot / 'site_local_quality.csv')
    sites = I.read_csv(snapshot / 'site_step_audit.csv')
    reviews, metadata = Q.upstream_reviews(sequences, candidates, mapped)
    expected = Q.aggregate_sites(prior, quality, reviews, metadata)
    I.require(sites == [{k: str(v) for k, v in r.items()} for r in expected], 'Corrected site audit inconsistent')
    I.require(len(quality) == b['mapped_options'] == 14310 and len(sites) == b['input_sites'] == 10799,
              'Context input accounting changed')
    key = lambda r: (r['site_id'], r['candidate_row'], r['label_asym_id'], r['model_id'])
    originals = {key(r): r for r in mapped if r['mapping_status'] == 'mapped'}
    I.require(len(originals) == len(quality) == len({key(r) for r in quality}), 'Mapped context option duplicated/lost')
    for r in quality:
        original = originals[key(r)]
        I.require(all(original[k] == r[k] for k in ('source', 'structure_file_id', 'structure_sha256', 'label_seq_id')),
                  'Context option/raw association changed')
    inputs += upstream + [Path(__file__).resolve(), Path(V.__file__).resolve()]
    return b, sequences, candidates, ledger, archive, quality, sites, originals, list(dict.fromkeys(inputs))


def atom_key(a):
    return (a.get('label_asym_id', ''), a.get('label_seq_id', ''), a.get('auth_asym_id', ''),
            a.get('auth_seq_id', ''), a.get('pdbx_PDB_ins_code', ''), a.get('label_comp_id', ''))


def atom_context(block, models):
    """All ASU/model atoms retained; no alternate or partner is discarded/prepared."""
    cells, residues, atom_ids = defaultdict(list), defaultdict(list), {}
    bad_coordinates = set()
    for a in M.records(block, '_atom_site.'):
        model = a.get('pdbx_PDB_model_num') or '1'
        if model not in models:
            continue
        xyz = tuple(Q.numeric(a.get(k)) for k in ('Cartn_x', 'Cartn_y', 'Cartn_z'))
        a = {k:a.get(k,'') for k in ('id','label_asym_id','label_seq_id','auth_asym_id','auth_seq_id',
                                    'pdbx_PDB_ins_code','label_comp_id','label_atom_id','label_alt_id',
                                    'type_symbol','occupancy')}
        a['model'] = model
        a['xyz'] = xyz
        residues[(model, atom_key(a))].append(a)
        identity = (model, a.get('id'))
        I.require(identity not in atom_ids, 'Duplicate atom ID within context model')
        atom_ids[identity] = a
        if any(v is None for v in xyz):
            bad_coordinates.add(model)
        elif a.get('type_symbol', '').upper() not in ('H', 'D') and Q.numeric(a.get('occupancy')) != 0:
            cells[(model,) + tuple(math.floor(v / 6.0) for v in xyz)].append(a)
    return cells, residues, atom_ids, bad_coordinates


def nearby(cells, model, point, target_key):
    center = tuple(math.floor(v / 6.0) for v in point)
    out = []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                for a in cells.get((model, center[0]+dx, center[1]+dy, center[2]+dz), []):
                    if atom_key(a) != target_key and sum((x-y)**2 for x, y in zip(a['xyz'], point)) <= 36:
                        out.append(a)
    return out


def coherent_atoms(atoms, alternate):
    # A blank SG cannot pick a nonblank side-chain conformer without additional evidence.
    selected = [a for a in atoms if not a.get('label_alt_id') or a.get('label_alt_id') == alternate]
    names = [a.get('label_atom_id') for a in selected if a.get('type_symbol', '').upper() not in ('H', 'D')]
    unresolved = (len(names) != len(set(names)) or any(Q.numeric(a.get('occupancy')) is None or
                  not 0 < Q.numeric(a.get('occupancy')) <= 1 or any(v is None for v in a['xyz']) for a in selected))
    return set(names), unresolved


def base_row(quality):
    return dict({k: quality[k] for k in ('site_id', 'protein_accession', 'canonical_cys_position', 'candidate_row',
                'source', 'structure_id', 'label_asym_id', 'model_id', 'structure_file_id', 'structure_sha256')},
                local_gate_status=quality['local_gate_status'], context_status='held', reason_code='',
                inherited_pae_status=quality.get('pae_status',''),
                inherited_max_target_neighbor_pae_bidirectional=quality.get('max_target_neighbor_pae_bidirectional',''),
                inherited_full_canonical_model=quality.get('full_canonical_model',''),
                inherited_fragment_boundary_review=quality.get('fragment_boundary_review',''),
                neighborhood_frame='predicted_model' if quality['source'] == 'AlphaFold_DB' else 'deposited_ASU',
                missing_coherent_target_atoms='', target_conformer_ambiguity='', neighbor_residue_count='',
                neighbor_other_subchains='', neighbor_nonpolymer_components='', neighbor_atom_gaps='',
                neighbor_alternate_or_occupancy_ambiguity='', neighbor_unrecognized_components='',
                declared_mutations_in_neighborhood='', entity_mutation_annotation='',
                candidate_assembly_ids='', assembly_definition_status='not_assessed',
                validation_status='not_applicable' if quality['source'] == 'AlphaFold_DB' else 'not_assessed',
                target_validation_attributes='', target_validation_children='',
                validation_neighbor_matched_count='', validation_neighbor_missing_count='',
                validation_report_sha256='', validation_revision_comparison='not_applicable',
                report_revision_date='', coordinate_revision_date='',
                structural_eligibility='held', chemical_state='not_assessed', protocol_version='0.1.0')


def evaluate_option(block, candidate, sequence, original, quality, atoms, tables, validation=None):
    row = base_row(quality)
    if quality['local_gate_status'] == 'excluded':
        row.update(context_status='not_assessed', reason_code='LOCAL_CONFIDENCE_OPTION_EXCLUDED')
        return row, []
    cells, residues, identities, bad = atoms
    model, label = original['model_id'], original['label_asym_id']
    sg = identities[(model, original['selected_SG_atom_id'])]
    I.require(sg['label_asym_id'] == label and sg['label_atom_id'] == 'SG' and
              sg['label_comp_id'] == 'CYS' and sg['type_symbol'].upper() == 'S' and
              all(abs(x-float(original[k])) < 1e-6 for x, k in zip(sg['xyz'], ('SG_x','SG_y','SG_z'))),
              'Frozen mapped SG/context identity mismatch')
    target = atom_key(sg)
    names, ambiguity = coherent_atoms(residues[(model, target)], original['selected_SG_alt_id'])
    near = nearby(cells, model, sg['xyz'], target)
    keys = sorted({atom_key(a) for a in near})
    gaps, unknown, alternates = [], [], []
    for key in keys:
        component = key[-1]
        group = residues[(model, key)]
        # Each deposited alternate is assessed separately; names are never unioned across alternates.
        choices = sorted({a.get('label_alt_id', '') for a in group if a.get('label_alt_id')}) or ['']
        if len(choices) > 1:
            alternates.append(list(key))
        if component not in SIDECHAINS:
            if component not in ('HOH', 'DOD'):
                unknown.append(list(key))
            continue
        expected = BACKBONE | set(SIDECHAINS[component].split())
        for alt in choices:
            found, unresolved = coherent_atoms(group, alt)
            if unresolved:
                alternates.append(list(key))
            if expected - found:
                gaps.append(dict(residue=list(key), alternate=alt, missing_atoms=sorted(expected-found)))
    entities, _, _, asym = M.polymer_index(block)
    entity = asym[label]
    mapping, _, changes = M.sequence_map(block, candidate, sequence, label, entity, entities[entity])
    reverse = {n: pos for pos, n in mapping.items()}
    neighbor_positions = {reverse[int(k[1])] for k in keys if k[0] == label and k[1].isdigit() and int(k[1]) in reverse}
    assemblies = sorted({g['assembly_id'] for g in tables['_pdbx_struct_assembly_gen.']
                         if label in [x.strip() for x in g.get('asym_id_list','').split(',')]})
    mutation = next((e.get('pdbx_mutation', '') for e in tables['_entity.'] if e.get('id') == entity), '')
    row.update(missing_coherent_target_atoms=';'.join(sorted(Q.CYS_ATOMS-names)), target_conformer_ambiguity=ambiguity,
               neighbor_residue_count=len(keys), neighbor_other_subchains=';'.join(sorted({k[0] for k in keys if k[0] != label})),
               neighbor_nonpolymer_components=';'.join(sorted({k[-1] for k in keys if not k[1]})),
               neighbor_atom_gaps=json.dumps(gaps, sort_keys=True),
               neighbor_unrecognized_components=json.dumps(unknown),
               neighbor_alternate_or_occupancy_ambiguity=json.dumps(alternates),
               declared_mutations_in_neighborhood=';'.join(str(n) for n in sorted(set(changes)&neighbor_positions)),
               entity_mutation_annotation=mutation, candidate_assembly_ids=';'.join(assemblies),
               assembly_definition_status='deposited_recipes_not_generated' if assemblies else 'no_matching_deposited_recipe',
               reason_code='STRUCTURAL_CONTEXT_REVIEW_PENDING' if model not in bad else 'MODEL_COORDINATE_CONTEXT_UNRESOLVED')
    revision_dates = [x.get('revision_date','') for x in tables['_pdbx_audit_revision_history.']]
    row['coordinate_revision_date'] = max(revision_dates, default='')
    evidence = []
    if quality['source'] == 'PDBe_SIFTS':
        if validation is None:
            row['validation_status'] = 'report_not_available_for_review'
        else:
            index, entry, report = validation
            status, payload = V.match_residue(index, model, sg.get('auth_asym_id',''), sg.get('auth_seq_id',''),
                                              sg.get('pdbx_PDB_ins_code',''), 'CYS', original['selected_SG_alt_id'])
            row.update(validation_status=status, validation_report_sha256=report['sha256'],
                       report_revision_date=entry.get('PDB-revision-date',''))
            left, right = row['coordinate_revision_date'], row['report_revision_date']
            row['validation_revision_comparison'] = 'dates_equal_manual_compatibility_pending' if left and left == right else 'dates_differ' if left and right else 'revision_comparison_unresolved'
            if payload is not None:
                row['target_validation_attributes'] = json.dumps(payload['attributes'], sort_keys=True)
                row['target_validation_children'] = json.dumps(payload['reported_children'], sort_keys=True)
            matched, missing = 0, 0
            for key in [target] + keys:
                group = residues[(model, key)]
                alts = sorted({a.get('label_alt_id','') for a in group}) or ['']
                for alt in alts:
                    state, details = V.match_residue(index, model, key[2], key[3], key[4], key[5], alt)
                    evidence.append(dict(structure_id=quality['structure_id'], report_sha256=report['sha256'], model_id=model,
                                         auth_chain=key[2], auth_seq_id=key[3], insertion_code=key[4], component=key[5], alternate=alt,
                                         validation_match_status=state, attributes=json.dumps(details['attributes'],sort_keys=True) if details else '',
                                         reported_children=json.dumps(details['reported_children'],sort_keys=True) if details else ''))
                    if key != target:
                        matched += state == 'matched_for_manual_review'
                        missing += state != 'matched_for_manual_review'
            row.update(validation_neighbor_matched_count=matched, validation_neighbor_missing_count=missing)
    return row, evidence


def acquire_reports(project, ledger, network, workers, timeout, attempts, probe=False):
    planned = V.plan(r['structure_id'] for r in ledger if r['kind'] == 'pdb_mmcif')
    plan_hash = V.sha(json.dumps(planned, sort_keys=True).encode())
    cache = project / 'data/raw/experimental_validation/step05b4' / plan_hash
    selected = planned[:1] if probe else planned
    def request(item):
        return V.download(cache, item, network=network, timeout=timeout, attempts=attempts)
    results = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for n, row in enumerate(pool.map(request, selected), 1):
            for field in ('local_file', 'metadata_file'):
                if row[field]:
                    row[field] = str(Path(row[field]).relative_to(project))
            results.append(row)
            if n % 25 == 0 or n == len(selected):
                print(f'  Validation reports checked: {n}/{len(selected)}', flush=True)
    return planned, results, plan_hash


def run(project, network=False, workers=2, timeout=120, attempts=4, probe=False):
    out = project / 'results/05b4_structure_context'
    if probe:
        out = out / 'probe'
    out.mkdir(parents=True, exist_ok=True)
    marker = 'PROBE_SUCCESS.txt' if probe else 'SUCCESS.txt'
    (out / marker).unlink(missing_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='context_stage_', dir=out))
    try:
        g = M.gemmi_module()
        b, seqs, candidates, ledger, archive, quality, sites, originals, inputs = load_input(project)
        print('[1/3] Acquiring/reusing official wwPDB experimental validation reports...', flush=True)
        planned, reports, plan_hash = acquire_reports(project, ledger, network or probe, workers, timeout, attempts, probe)
        I.write_csv(out/'validation_download_plan.csv', planned, list(planned[0]))
        I.write_csv(out/'validation_download_manifest.csv', reports, list(reports[0]))
        failures = [r for r in reports if r['status'] == 'failed']
        (out/'failure_context.json').write_text(json.dumps(failures,indent=2)+'\n')
        I.require(not failures, 'Validation acquisition failed; see validation_download_manifest.csv and failure_context.json')
        if probe:
            I.require(reports[0]['status'] == 'downloaded', 'Probe report not successfully acquired')
            (out/'FAILURE.txt').unlink(missing_ok=True)
            (out/marker).write_text('PROBE_VALIDATION_COMPLETE\n')
            print('PROBE_VALIDATION_COMPLETE: ' + reports[0]['structure_id'], flush=True)
            return dict(status='PROBE_VALIDATION_COMPLETE')
        print('[2/3] Verifying all original raw files and collecting context evidence...', flush=True)
        verified = M.verify_raw(project, ledger, archive)
        key = lambda r: (r['site_id'], r['candidate_row'], r['label_asym_id'], r['model_id'])
        groups = defaultdict(list)
        for r in quality:
            groups[r['structure_file_id']].append(r)
        report_by_entry = {r['structure_id']: r for r in reports}
        rows, entries, annotations, evidence, validation_entries = [], [], [], [], []
        coordinates = [r for r in ledger if r['kind'].endswith('mmcif')]
        for number, raw in enumerate(coordinates, 1):
            block = g.cif.read_file(str(I.safe_path(project, raw['local_file']))).sole_block()
            I.require(block.name.casefold() == raw['structure_id'].casefold(), 'Context coordinate identity mismatch')
            tables = {category:list(M.records(block,category)) for category in ANNOTATIONS}
            for category, data in tables.items():
                for ordinal, record in enumerate(data, 1):
                    annotations.append(dict(structure_id=raw['structure_id'], structure_file_id=raw['file_id'],
                                            structure_sha256=raw['sha256'], category=category, row_number=ordinal,
                                            deposited_record=json.dumps(record,sort_keys=True)))
            entries.append(dict(structure_id=raw['structure_id'],source_kind=raw['kind'],structure_file_id=raw['file_id'],
                                structure_sha256=raw['sha256'],assembly_definitions=len(tables['_pdbx_struct_assembly.']),
                                assembly_generators=len(tables['_pdbx_struct_assembly_gen.']),
                                assembly_operators=len(tables['_pdbx_struct_oper_list.']),
                                deposited_connections=len(tables['_struct_conn.']),
                                declared_reference_differences=len(tables['_struct_ref_seq_dif.']),
                                reported_unobserved_residues=len(tables['_pdbx_unobs_or_zero_occ_residues.']),
                                reported_unobserved_atoms=len(tables['_pdbx_unobs_or_zero_occ_atoms.']),
                                final_structural_pass=False,protocol_version='0.1.0'))
            options = groups[raw['file_id']]
            assessed = [r for r in options if r['local_gate_status'] != 'excluded']
            atoms = atom_context(block,{r['model_id'] for r in assessed}) if assessed else None
            report = report_by_entry.get(raw['structure_id'])
            validation = None
            if report is not None and report['status'] == 'downloaded':
                body = I.safe_path(project,report['local_file'])
                I.require(M.file_sha(body) == report['sha256'], 'Validation raw hash changed after acquisition')
                root, entry = V.xml_content(body.read_bytes(), raw['structure_id'])
                validation = V.residue_index(root), entry, report
                validation_entries.append(dict(structure_id=raw['structure_id'],report_sha256=report['sha256'],
                                                schema_attributes=json.dumps(dict(root.attrib),sort_keys=True),
                                                entry_attributes=json.dumps(entry,sort_keys=True)))
                inputs += [body, I.safe_path(project,report['metadata_file'])]
            for q in options:
                original = originals[key(q)]
                candidate = candidates[int(q['candidate_row'])-1]
                result, used = evaluate_option(block,candidate,seqs[q['protein_accession']],original,q,atoms,tables,validation)
                rows.append(result)
                evidence += used
            del atoms, block, validation
            if number % 25 == 0 or number == len(coordinates):
                print(f'  Coordinate contexts checked: {number}/{len(coordinates)}',flush=True)
        I.require(len(rows) == len(quality), 'Context mapped-option accounting lost')
        rows.sort(key=key)
        # Preserve inherited terminal uncertainty exactly; this evidence inventory issues no context pass.
        audit = [dict(r, step_id='05b4', reason_code='CONTEXT_REVIEW_PENDING' if r['status']=='held' else r['reason_code']) for r in sites]
        queue = {}
        for r in rows:
            if r['local_gate_status'] == 'excluded':
                continue
            k = (r['candidate_row'],r['label_asym_id'],r['model_id'])
            if k not in queue:
                queue[k] = dict(candidate_row=r['candidate_row'],protein_accession=r['protein_accession'],source=r['source'],
                                structure_id=r['structure_id'],label_asym_id=r['label_asym_id'],model_id=r['model_id'],
                                assembly_options=r['candidate_assembly_ids'],review_status='pending',reviewer='',evidence_note='',
                                sites_requiring_review=0,structural_eligibility='held',protocol_version='0.1.0')
            queue[k]['sites_requiring_review'] += 1
        deduplicated = {json.dumps(r,sort_keys=True):r for r in evidence}
        outputs = [('site_context_evidence.csv',rows),('structure_context_inventory.csv',entries),
                   ('context_annotations.csv',annotations),('site_step_audit.csv',audit),
                   ('manual_review_queue.csv',list(queue.values())),('raw_file_checksums.csv',verified)]
        for name, data in outputs:
            I.require(bool(data),'Unexpected empty context output')
            I.write_csv(stage/name,data,list(data[0]))
        evidence_fields = ['structure_id','report_sha256','model_id','auth_chain','auth_seq_id','insertion_code','component',
                           'alternate','validation_match_status','attributes','reported_children']
        I.write_csv(stage/'experimental_validation_residues.csv',list(deduplicated.values()),evidence_fields)
        I.write_csv(stage/'validation_entry_inventory.csv',validation_entries,
                    ['structure_id','report_sha256','schema_attributes','entry_attributes'])
        inputs = list(dict.fromkeys(inputs))
        I.write_csv(stage/'input_checksums.csv',[dict(input=str(f.relative_to(project)),sha256=M.file_sha(f)) for f in inputs],
                    ['input','sha256'])
        summary = dict(status='STEP05B4_GENERATED_REVIEW_PENDING',implementation_version=VERSION,protocol_version='0.1.0',
                       input_proteins=len(seqs),input_sites=len(audit),mapped_options=len(rows),coordinate_files=len(entries),
                       context_options_assessed=sum(r['local_gate_status']!='excluded' for r in rows),
                       validation_plan_sha256=plan_hash,validation_report_status_counts=dict(Counter(r['status'] for r in reports)),
                       experimental_target_validation_status_counts=dict(Counter(r['validation_status'] for r in rows if r['source']=='PDBe_SIFTS')),
                       inherited_site_status_counts=dict(Counter(r['status'] for r in audit)),review_queue_rows=len(queue),
                       final_structural_passes=0,python=platform.python_version(),gemmi=g.__version__,
                       execution_arguments=dict(download_validation=network,workers=workers,timeout=timeout,attempts=attempts),
                       generated_utc=datetime.now(timezone.utc).isoformat())
        (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        report = '\n'.join(['# Step 05b4 — structural context evidence inventory','', 'Status: '+summary['status'],'',
                           f'- Input proteins: {len(seqs)}; canonical Cys sites: {len(audit)}',
                           f'- Coordinate files: {len(entries)}; mapped options retained: {len(rows)}',
                           '- Official validation report statuses: '+str(summary['validation_report_status_counts']),
                           '- Experimental target report matches: '+str(summary['experimental_target_validation_status_counts']),
                           '- Inherited site statuses: '+str(summary['inherited_site_status_counts']),
                           f'- Review queue rows: {len(queue)}; final structural passes: 0','',
                           'Experimental neighborhoods refer to deposited ASUs, not generated/selected biological assemblies.',
                           'All deposited assembly recipes/operators, connections, sequence differences and unobserved-atom/residue annotations remain evidence.',
                           'Target coherent atoms and observed modeled neighbors are audited without atom repair or alternate-state merging.',
                           'Official XML validation fields/children remain raw: missing metrics are not zero and no new resolution/RSRZ cutoff is imposed.',
                           'Report/coordinate revision dates are compared; equality alone does not prove identical validated coordinates.',
                           'Native, interdomain, assembly and chemical context still require review; no structure is selected for features.',
                           'No SASA, pKa or oxidation score is generated.',''])
        (stage/'context_report.md').write_text(report)
        for f in stage.iterdir():
            f.replace(out/f.name)
        (out/'FAILURE.txt').unlink(missing_ok=True)
        (out/marker).write_text(summary['status']+'\n')
        print('[3/3] Outputs: '+str(out),flush=True)
        print(report)
        return summary
    except Exception as e:
        (out/'FAILURE.txt').write_text(str(e)+'\n')
        raise
    finally:
        shutil.rmtree(stage,ignore_errors=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project',type=Path,default=Path(__file__).resolve().parents[1])
    parser.add_argument('--download-validation',action='store_true',help='Download uncached official XML reports; otherwise use cache only')
    parser.add_argument('--probe-validation',action='store_true',help='Test one official report before full execution')
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--timeout',type=int,default=120)
    parser.add_argument('--attempts',type=int,default=4)
    args = parser.parse_args()
    if not (1<=args.workers<=4 and 1<=args.timeout<=600 and 1<=args.attempts<=8):
        parser.error('Invalid workers/timeout/attempts')
    try:
        run(args.project.resolve(),args.download_validation,args.workers,args.timeout,args.attempts,args.probe_validation)
    except Exception as e:
        print('STEP05B4 FAILED: '+str(e),file=sys.stderr)
        sys.exit(1)
