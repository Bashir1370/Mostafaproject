"""Structural availability is not local eligibility; identity and errors matter."""
import importlib.util
import io
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch
import urllib.error

PROJECT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('catalogue', PROJECT / 'scripts/05a_catalogue_GSE286387.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)


def model(**changes):
    result = dict(modelEntityId='AF-P12345-F1', sequence='MCC', sequenceStart=1, sequenceEnd=3,
                  latestVersion=6, taxId=10090, uniprotAccession='P12345', chainId='A',
                  cifUrl='https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v6.cif')
    result.update(changes)
    return result


class CatalogueIdentity(unittest.TestCase):
    def test_current_schema_and_full_or_partial_sequence(self):
        r = C.parse_afdb('P12345', 'MCC', [model()])[0]
        self.assertEqual(r['sequence_match'], 'exact_full')
        self.assertEqual(r['candidate_status'], 'candidate')
        r = C.parse_afdb('P12345', 'MCCAAA', [model(sequence='CC', sequenceStart=2, sequenceEnd=3)])[0]
        self.assertEqual(r['sequence_match'], 'exact_partial')
        self.assertIn('BOUNDARY_REVIEW_PENDING', r['reason'])

    def test_isoform_taxonomy_and_sequence_mismatch_are_held(self):
        for changes in [dict(uniprotAccession='P12345-2'), dict(taxId=9606), dict(sequence='MCA'),
                        dict(sequenceEnd=4), dict(isComplex=True)]:
            self.assertEqual(C.parse_afdb('P12345', 'MCC', [model(**changes)])[0]['candidate_status'], 'held')

    def test_conflicting_schema_or_missing_download_url_not_accepted(self):
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            C.parse_afdb('P12345', 'MCC', [model(uniprotSequence='AAA')])
        self.assertEqual(C.parse_afdb('P12345', 'MCC', [model(cifUrl='')])[0]['candidate_status'], 'held')

    def test_pdbe_retains_chains_and_taxonomy_is_not_global_resolution(self):
        mouse = dict(pdb_id='1abc', chain_id='A', unp_start=1, unp_end=3, tax_id=10090,
                     resolution=9.9, experimental_method='X-ray diffraction', coverage=1)
        human = dict(mouse, chain_id='B', tax_id=9606, resolution=1.0)
        unknown = dict(mouse, chain_id='C', tax_id=None)
        long_interval = dict(mouse, chain_id='D', unp_end=4)
        rows = C.parse_pdbe('P12345', {'P12345': [mouse, human, unknown, long_interval, mouse]}, 3)
        self.assertEqual(len(rows), 4)
        self.assertEqual([r['candidate_status'] for r in rows], ['candidate', 'excluded', 'held', 'held'])
        self.assertEqual(rows[0]['sequence_match'], 'not_assessed')

    def test_frozen_inventory_is_loaded_before_network(self):
        binding, sequences, inputs = C.load_input(PROJECT)
        self.assertEqual(len(sequences), 760)
        self.assertEqual(sum(s.count('C') for s in sequences.values()), 10799)

    def test_changed_frozen_input_fails_before_network_and_removes_old_success(self):
        binding_path = PROJECT / 'config/GSE286387_step05_input.json'
        binding = json.loads(binding_path.read_text())
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            for relative in list(binding['input_sha256']) + ['config/GSE286387_step05_input.json']:
                target = project / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(PROJECT / relative, target)
            fasta = project / binding['fasta_file']
            fasta.write_bytes(fasta.read_bytes() + b'\n')
            out = project / 'results/05a_structure_catalogue'
            out.mkdir(parents=True)
            (out / 'SUCCESS.txt').write_text('old success')
            with patch.object(C, 'fetch', side_effect=AssertionError('network called')):
                with self.assertRaisesRegex(ValueError, 'Frozen Step 05 input changed'):
                    C.run(project)
            self.assertFalse((out / 'SUCCESS.txt').exists())
            self.assertTrue((out / 'FAILURE.txt').exists())

    def test_unresolved_offered_model_is_held_not_a_missing_record(self):
        provenance = dict(url='test', sha256='test', http_status=200, retrieved_utc='test', service_release='test')
        with patch.object(C, 'fetch', side_effect=[(None, provenance), ([model(uniprotAccession='P12345-2')], provenance)]):
            self.assertEqual(C.catalogue_one('P12345', 'MCC', Path('.'), False)[1]['status'], 'held')
        with patch.object(C, 'fetch', side_effect=[(None, provenance), (None, provenance)]):
            self.assertEqual(C.catalogue_one('P12345', 'MCC', Path('.'), False)[1]['status'], 'unassessable')


class CatalogueCache(unittest.TestCase):
    def test_forbidden_names_endpoint_and_is_not_cached_or_retried(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp); url = C.AFDB + 'P12345'
            error = urllib.error.HTTPError(url, 403, 'Forbidden',
                                          {'Server': 'example', 'Set-Cookie': 'must-not-log'},
                                          io.BytesIO(b'access denied'))
            with patch.object(C.urllib.request, 'urlopen', side_effect=error) as request:
                with self.assertRaises(C.CatalogueHTTPError) as captured:
                    C.fetch(url, cache)
            self.assertEqual(request.call_count, 1)
            self.assertIn(url, str(captured.exception))
            self.assertEqual(captured.exception.details['response_excerpt'], 'access denied')
            self.assertNotIn('Set-Cookie', captured.exception.details['response_headers'])
            self.assertFalse(list(cache.iterdir()))

    def test_http_failure_writes_diagnostic_and_invalidates_success(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            original = urllib.error.HTTPError(C.PDBE + 'P12345', 403, 'Forbidden', {}, io.BytesIO(b'denied'))
            error = C.CatalogueHTTPError(C.PDBE + 'P12345', original)
            out = project / 'results/05a_structure_catalogue'
            out.mkdir(parents=True)
            (out / 'SUCCESS.txt').write_text('old success')
            with patch.object(C, 'load_input', side_effect=error):
                with self.assertRaises(C.CatalogueHTTPError):
                    C.run(project)
            self.assertFalse((out / 'SUCCESS.txt').exists())
            self.assertTrue((out / 'FAILURE.txt').exists())
            diagnostic = json.loads((out / 'failure_context.json').read_text())
            self.assertEqual(diagnostic['http_status'], 403)
            self.assertEqual(diagnostic['requested_url'], C.PDBE + 'P12345')

    def test_cached_404_is_not_a_network_error_and_corruption_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp); url = C.AFDB + 'P12345'; key = C.I.sha(url.encode()); raw = b'not found'
            (cache / (key + '.response')).write_bytes(raw)
            meta = dict(url=url, sha256=C.I.sha(raw), http_status=404, retrieved_utc='test', service_release='not_supplied_by_service')
            (cache / (key + '.metadata.json')).write_text(json.dumps(meta))
            with patch.object(C.urllib.request, 'urlopen', side_effect=AssertionError('network used')):
                self.assertIsNone(C.fetch(url, cache, True)[0])
            (cache / (key + '.response')).write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'Corrupt'):
                C.fetch(url, cache, True)

    def test_server_error_never_becomes_missing_structure(self):
        with tempfile.TemporaryDirectory() as temp:
            cache = Path(temp); url = C.AFDB + 'P12345'
            error = urllib.error.HTTPError(url, 503, 'unavailable', {}, io.BytesIO(b'error'))
            with patch.object(C.urllib.request, 'urlopen', side_effect=error), patch.object(C.time, 'sleep'):
                with self.assertRaises(urllib.error.HTTPError):
                    C.fetch(url, cache)
            self.assertFalse(list(cache.iterdir()))

    def test_offline_missing_cache_is_explicit_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaisesRegex(ValueError, 'Offline cache missing'):
                C.fetch(C.AFDB + 'P12345', Path(temp), True)


@unittest.skipUnless((PROJECT / 'results/05a_structure_catalogue/SUCCESS.txt').exists(), 'Run live catalogue first')
class RealCatalogue(unittest.TestCase):
    def test_every_protein_and_both_queries_are_accounted_for_without_site_approval(self):
        out = PROJECT / 'results/05a_structure_catalogue'
        proteins = C.I.read_csv(out / 'protein_structure_inventory.csv')
        sources = C.I.read_csv(out / 'source_manifest.csv')
        candidates = C.I.read_csv(out / 'structure_candidates.csv')
        self.assertEqual(len(proteins), 760)
        self.assertEqual(len({r['protein_accession'] for r in proteins}), 760)
        self.assertEqual(len(sources), 1520)
        self.assertEqual(len({r['url'] for r in sources}), 1520)
        self.assertTrue(all(r['structural_eligibility'] == 'not_assessed' for r in candidates))
        self.assertTrue(all(r['status'] in ('held', 'unassessable') for r in proteins))


class FrozenReference(unittest.TestCase):
    def test_reference_run_without_any_network_and_same_audits(self):
        binding, sequences, _ = C.load_input(PROJECT)
        cache, ref, inputs = C.reference_cache(PROJECT, binding, sequences)
        self.assertEqual(len(list(cache.iterdir())), 3040)
        with patch.object(C.urllib.request, 'urlopen', side_effect=AssertionError('Network forbidden')):
            summary = C.run(PROJECT, workers=4, use_reference_cache=True)
        self.assertEqual(summary['metadata_mode'], 'frozen_reference')
        self.assertTrue(summary['offline'])
        self.assertEqual(summary['query_records'], 1520)
        self.assertEqual(summary['candidate_status_counts'], {'held': 324, 'candidate': 1533, 'excluded': 15})
        self.assertEqual(C.I.sha((PROJECT / 'results/05a_structure_catalogue/source_manifest.csv').read_bytes()),
                         ref['source_manifest_sha256'])

    def test_bad_reference_hash_fails_before_import(self):
        binding, sequences, _ = C.load_input(PROJECT)
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            for relative in ['references/step05a_metadata_2026-10-09',
                             'docs/audits/Structure_catalogue_validation']:
                shutil.copytree(PROJECT / relative, project / relative)
            descriptor = project / 'references/step05a_metadata_2026-10-09/reference.json'
            ref = json.loads(descriptor.read_text()); ref['archive_sha256'] = '0' * 64
            descriptor.write_text(json.dumps(ref))
            with self.assertRaisesRegex(ValueError, 'Reference checksum'):
                C.reference_cache(project, binding, sequences)
            self.assertFalse((project / 'data').exists())

    def test_unsafe_archive_is_rejected_even_with_updated_checksum(self):
        binding, sequences, _ = C.load_input(PROJECT)
        import tarfile
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            for relative in ['references/step05a_metadata_2026-10-09',
                             'docs/audits/Structure_catalogue_validation']:
                shutil.copytree(PROJECT / relative, project / relative)
            descriptor = project / 'references/step05a_metadata_2026-10-09/reference.json'
            ref = json.loads(descriptor.read_text())
            archive = project / ref['archive_file']
            with tarfile.open(archive, 'w:gz') as bundle:
                member = tarfile.TarInfo('../escape'); member.size = 1
                bundle.addfile(member, io.BytesIO(b'x'))
            ref['archive_sha256'] = C.I.sha(archive.read_bytes()); descriptor.write_text(json.dumps(ref))
            with self.assertRaisesRegex(ValueError, 'Unsafe/unexpected/duplicate'):
                C.reference_cache(project, binding, sequences)
            self.assertFalse((project / 'data').exists())

    def test_existing_corrupt_reference_cache_is_not_reused(self):
        binding, sequences, _ = C.load_input(PROJECT)
        cache, ref, _ = C.reference_cache(PROJECT, binding, sequences)
        target = next(cache.glob('*.response')); original = target.read_bytes()
        try:
            target.write_bytes(b'changed')
            with self.assertRaisesRegex(ValueError, 'Existing reference cache changed'):
                C.reference_cache(PROJECT, binding, sequences)
        finally:
            target.write_bytes(original)


if __name__ == '__main__':
    unittest.main()
