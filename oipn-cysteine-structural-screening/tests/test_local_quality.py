"""Confidence boundary, source consistency, geometry and complete audit contracts."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

PROJECT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('quality', PROJECT / 'scripts/05b3_local_quality_GSE286387.py')
Q = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(Q)
OTHER = importlib.util.spec_from_file_location('mapping_test', Path(__file__).with_name('test_structure_mapping.py'))
T = importlib.util.module_from_spec(OTHER)
OTHER.loader.exec_module(T)
FIXED = dict(local_radius_angstrom=6.0, predicted_target_plddt_minimum=90,
             predicted_neighbor_plddt_minimum=70, predicted_neighbor_plddt_strict=90)


def model(target=90, neighbor=70, distant=10, distance=6, omit='', alt='.', occupancy=1,
          extra_element='C', sequence='ACG', start=1):
    lines = []
    for k, (name, element) in enumerate([('N', 'N'), ('CA', 'C'), ('C', 'C'), ('O', 'O'), ('CB', 'C'), ('SG', 'S')], 1):
        if name == omit:
            continue
        lines.append(f'{k} {element} {name} {alt} CYS L 2 A 52 B 1 {occupancy} 0 0 0 {target}\n')
    lines.append(f'7 C CA . ALA L 1 A 51 . 1 1 {distance} 0 0 {neighbor}\n')
    lines.append(f'8 {extra_element} CA . GLY L 3 A 53 . 1 1 6.001 0 0 {distant}\n')
    block = T.fixture(atoms=''.join(lines))
    candidate = T.candidate('AlphaFold_DB')
    candidate.update(sequence_start=str(start), sequence_end=str(start + 2))
    mapped = Q.M.map_candidate(block, candidate, 1, sequence, T.FILE)
    conf = dict(residueNumber=[1, 2, 3], confidenceScore=[neighbor, target, distant])
    pae = [dict(predicted_aligned_error=[[0, 2, 3], [1, 0, 4], [3, 2, 0]], max_predicted_aligned_error=31.75)]
    return block, candidate, mapped, conf, pae


def evaluate(*args, **kwargs):
    block, candidate, mapped, conf, pae = model(*args, **kwargs)
    sequence = kwargs.get('sequence', 'ACG')
    return Q.evaluate_model(block, candidate, sequence, mapped, conf, pae, FIXED)[0]


class LocalQuality(unittest.TestCase):
    def test_inclusive_thresholds_and_radius_ignore_distant_low_confidence(self):
        row = evaluate()
        self.assertEqual(row['local_gate_status'], 'pass')
        self.assertEqual(row['neighbor_residue_count'], 1)
        self.assertEqual(row['minimum_neighbor_plddt'], 70)
        self.assertEqual(row['max_target_neighbor_pae_bidirectional'], 2)
        self.assertEqual(row['structural_eligibility'], 'held')

    def test_below_target_and_neighbor_thresholds(self):
        self.assertEqual(evaluate(target=89.99)['reason_code'], 'TARGET_PLDDT_BELOW_90')
        self.assertEqual(evaluate(neighbor=69.99)['reason_code'], 'NEIGHBOR_PLDDT_BELOW_70')

    def test_radius_excludes_just_outside_and_target_atoms(self):
        row = evaluate(neighbor=10, distance=6.001)
        self.assertEqual(row['local_gate_status'], 'pass')
        self.assertEqual(row['neighbor_residue_count'], 0)
        self.assertEqual(row['minimum_neighbor_plddt'], '')

    def test_hydrogen_and_deuterium_are_not_neighbor_atoms(self):
        for element in ('H', 'D'):
            block, candidate, mapped, conf, pae = model(extra_element=element)
            table = block.find_mmcif_category('_atom_site.')
            # Move the low-confidence H/D atom inside the radius.
            table[-1][list(table.tags).index('_atom_site.Cartn_x')] = '1'
            row = Q.evaluate_model(block, candidate, 'ACG', mapped, conf, pae, FIXED)[0]
            self.assertEqual(row['minimum_neighbor_plddt'], 70)

    def test_coherent_target_requires_all_heavy_atoms(self):
        self.assertEqual(evaluate(omit='CB')['reason_code'], 'COHERENT_TARGET_HEAVY_ATOMS_MISSING')

    def test_sensitivity_does_not_change_primary_gate(self):
        self.assertFalse(evaluate()['neighbor_90_sensitivity'])
        self.assertTrue(evaluate(neighbor=90)['neighbor_90_sensitivity'])

    def test_cif_json_disagreement_holds(self):
        block, c, m, confidence, pae = model()
        confidence['confidenceScore'][1] = 99
        row = Q.evaluate_model(block, c, 'ACG', m, confidence, pae, FIXED)[0]
        self.assertEqual(row['reason_code'], 'CIF_JSON_CONFIDENCE_DISAGREEMENT')

    def test_confidence_numbering_nan_and_boolean_hold(self):
        for data in [dict(residueNumber=[2, 1, 3], confidenceScore=[90]*3),
                     dict(residueNumber=[1, 2, 3], confidenceScore=[90, float('nan'), 90]),
                     dict(residueNumber=[1, 2, 3], confidenceScore=[90, True, 90])]:
            with self.assertRaises(Q.QualityHold):
                Q.confidence_values(data, 3)

    def test_alternates_and_nonunit_occupancy_hold(self):
        self.assertEqual(evaluate(alt='A')['reason_code'], 'PREDICTED_ALTERNATE_OR_OCCUPANCY_CONTEXT_UNRESOLVED')
        self.assertEqual(evaluate(occupancy=.5)['local_gate_status'], 'held')

    def test_fragment_uses_model_relative_confidence_and_retains_boundary_flag(self):
        row = evaluate(sequence='MMACGMM', start=3)
        self.assertEqual(row['local_gate_status'], 'pass')
        self.assertEqual(row['canonical_cys_position'], 4)
        self.assertEqual(row['neighbor_canonical_positions'], '3')
        self.assertTrue(row['fragment_boundary_review'])

    def test_pae_rounding_and_missing_matrix_not_imputed_or_new_gate(self):
        matrix = Q.pae_matrix([dict(predicted_aligned_error=[[0, 32], [32, 0]], max_predicted_aligned_error=31.75)], 2)
        self.assertEqual(matrix[0][1], 32)
        block, c, m, conf, _ = model()
        row = Q.evaluate_model(block, c, 'ACG', m, conf, [], FIXED)[0]
        self.assertEqual(row['local_gate_status'], 'pass')
        self.assertEqual(row['pae_status'], 'PAE_SCHEMA_UNRESOLVED')
        self.assertEqual(row['max_target_neighbor_pae_bidirectional'], '')

    def test_mapped_coordinate_disagreement_fails_instead_of_scientific_absence(self):
        block, c, m, conf, pae = model()
        m[0]['SG_x'] = 99
        with self.assertRaisesRegex(ValueError, 'MAPPED_SG_COORDINATES_DISAGREE'):
            Q.evaluate_model(block, c, 'ACG', m, conf, pae, FIXED)

    def test_aggregation_preserves_unresolved_and_no_structure_sites(self):
        audit = Q.M.site_audit({'P': 'CC', 'Q': 'C'}, [], {'P'})
        rows = [dict(site_id='P:C1', local_gate_status='pass')]
        result = Q.aggregate_sites(audit, rows)
        self.assertEqual([r['status'] for r in result], ['held', 'held', 'unassessable'])
        self.assertEqual(result[1]['structural_eligibility'], 'held')

    def test_full_frozen_mapping_binding_and_all_raw_hash_gate(self):
        values = Q.load_input(PROJECT)
        self.assertEqual(len(values[6]), 25026)
        with tempfile.TemporaryDirectory() as t:
            with self.assertRaisesRegex(ValueError, 'missing or changed'):
                Q.M.verify_raw(Path(t), values[4], values[5])

    def test_real_sample_agrees_with_independent_brute_force_neighbors(self):
        _, seqs, candidates, _, ledger, _, mapped, _, _, _ = Q.load_input(PROJECT)
        f = {r['kind']: r for r in ledger if r['structure_id'] == 'AF-A0A087WRH0-F1'}
        if any(not (PROJECT / r['local_file']).is_file() for r in f.values()):
            self.skipTest('Representative cached raw bytes not available')
        for r in f.values():
            self.assertEqual(Q.M.file_sha(PROJECT / r['local_file']), r['sha256'])
        block = Q.M.gemmi_module().cif.read_file(str(PROJECT / f['alphafold_mmcif']['local_file'])).sole_block()
        options = [r for r in mapped if r['structure_id'] == 'AF-A0A087WRH0-F1' and r['mapping_status'] == 'mapped']
        c = candidates[int(options[0]['candidate_row']) - 1]
        scores = json.loads((PROJECT / f['alphafold_confidence']['local_file']).read_text())
        pae = json.loads((PROJECT / f['alphafold_pae']['local_file']).read_text())
        result = Q.evaluate_model(block, c, seqs[c['protein_accession']], options, scores, pae, FIXED)
        atoms = list(Q.M.records(block, '_atom_site.'))
        confidence = dict(zip(scores['residueNumber'], scores['confidenceScore']))
        for original, row in zip(options, result):
            n = int(original['label_seq_id'])
            point = [float(original[k]) for k in ('SG_x', 'SG_y', 'SG_z')]
            near = {int(a['label_seq_id']) for a in atoms if a['type_symbol'].upper() not in ('H', 'D') and
                    int(a['label_seq_id']) != n and sum((float(a[k]) - x)**2 for k, x in
                    zip(('Cartn_x', 'Cartn_y', 'Cartn_z'), point)) <= 36}
            expected = confidence[n] >= 90 and all(confidence[i] >= 70 for i in near)
            self.assertEqual(row['local_gate_status'] == 'pass', expected)
            self.assertEqual(row['neighbor_residue_count'], len(near))
            self.assertEqual(row['pae_status'], 'parsed_descriptive_only')
        self.assertEqual(sum(r['local_gate_status'] == 'pass' for r in result), 21)

    def test_end_to_end_preserves_all_sites_and_experimental_options(self):
        b, seqs, candidates, links, ledger, archive, mapped, _, fixed, _ = Q.load_input(PROJECT)
        selected = [r for r in mapped if r['mapping_status'] == 'mapped' and
                    r['structure_id'] in ('AF-A0A087WRH0-F1', '11gl')]
        small_ledger = [r for r in ledger if r['structure_id'] in ('AF-A0A087WRH0-F1', '11gl')]
        if any(not (PROJECT / r['local_file']).is_file() for r in small_ledger):
            self.skipTest('Representative cached raw bytes not available')
        audit = Q.M.site_audit(seqs, selected, {c['protein_accession'] for c in candidates if c['candidate_status'] == 'held'})
        audit = [{k: str(v) for k, v in r.items()} for r in audit]
        with tempfile.TemporaryDirectory() as t:
            project = Path(t)
            for r in small_ledger:
                body = project / r['local_file']
                body.parent.mkdir(parents=True, exist_ok=True)
                body.write_bytes((PROJECT / r['local_file']).read_bytes())
                meta = dict(file_id=r['file_id'], url=r['url'], response_url=r['response_url'], http_status=200,
                            bytes=int(r['bytes']), sha256=r['sha256'], retrieved_utc=r['retrieved_utc'],
                            etag=r['etag'], last_modified=r['last_modified'])
                (project / r['metadata_file']).write_text(json.dumps(meta, indent=2) + '\n')
            values = b, seqs, candidates, links, small_ledger, archive, selected, audit, fixed, []
            with patch.object(Q, 'load_input', return_value=values):
                summary = Q.run(project)
            out = project / 'results/05b3_local_quality'
            self.assertEqual(summary['input_sites'], 10799)
            self.assertEqual(summary['mapped_options'], 43)
            self.assertEqual(summary['sites_with_local_confidence_options'], 21)
            self.assertEqual(summary['local_gate_status_counts'], {'held': 8, 'pass': 21, 'excluded': 14})
            self.assertEqual(len(Q.I.read_csv(out / 'site_step_audit.csv')), 10799)
            self.assertTrue((out / 'SUCCESS.txt').exists())
            self.assertEqual(summary['final_structural_passes'], 0)
            # A raw-byte mismatch invalidates SUCCESS; old tables cannot count as a complete run.
            (project / small_ledger[0]['local_file']).write_bytes(b'changed')
            with patch.object(Q, 'load_input', return_value=values):
                with self.assertRaises(ValueError):
                    Q.run(project)
            self.assertFalse((out / 'SUCCESS.txt').exists())
            self.assertTrue((out / 'FAILURE.txt').exists())


if __name__ == '__main__':
    unittest.main()
