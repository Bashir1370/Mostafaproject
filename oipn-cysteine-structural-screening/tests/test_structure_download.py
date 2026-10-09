"""Archive/transport contracts do not imply structural or site eligibility."""
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import urllib.error

PROJECT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('download', PROJECT / 'scripts/05b1_download_GSE286387.py')
D = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(D)
RAW = b'data_1abc\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n#\n'


def item():
    url = 'https://files.rcsb.org/download/1abc.cif'
    return dict(file_id=D.I.sha(('pdb_mmcif\n'+url).encode()), kind='pdb_mmcif',
                url=url, structure_id='1abc', required=True)


class Response(io.BytesIO):
    status = 200
    def __init__(self, raw=RAW, url=None):
        super().__init__(raw); self.headers={'Content-Length':str(len(raw))}; self.url=url or item()['url']
    def geturl(self): return self.url


class DownloadContracts(unittest.TestCase):
    def test_real_frozen_input_plan_and_complete_candidate_accounting(self):
        binding, rows, inputs = D.load_input(PROJECT)
        plan, links = D.make_plan(rows)
        self.assertEqual(len(plan), 2638)
        self.assertEqual(sum(r['required'] for r in plan), 1176)
        self.assertEqual(len(links), 1872)
        self.assertEqual(len({r['file_id'] for r in plan}),2638)
        self.assertEqual(sum(r['acquisition_action']=='planned' for r in links),1533)
        self.assertTrue(all(r['structural_eligibility']=='not_assessed' for r in links))

    def test_atomic_download_and_offline_reuse_without_network(self):
        with tempfile.TemporaryDirectory() as t:
            cache=Path(t)
            with patch.object(D.urllib.request,'urlopen',return_value=Response()):r=D.download(cache,item())
            self.assertEqual(r['status'],'downloaded')
            with patch.object(D.urllib.request,'urlopen',side_effect=AssertionError('network')):
                replay=D.download(cache,item(),offline=True)
            self.assertEqual(replay['sha256'],D.I.sha(RAW))
            self.assertEqual(replay['retrieved_utc'],r['retrieved_utc'])
            body,_=D.paths(cache,item());body.write_bytes(b'changed')
            self.assertEqual(D.download(cache,item(),offline=True)['status'],'failed')

    def test_html_wrong_identifier_and_truncated_payload_never_cached(self):
        for raw in (b'<html>access denied</html>', RAW.replace(b'1abc',b'2abc')):
            with tempfile.TemporaryDirectory() as t, patch.object(D.urllib.request,'urlopen',return_value=Response(raw)):
                cache=Path(t);self.assertEqual(D.download(cache,item())['status'],'failed')
                self.assertFalse(list(cache.iterdir()))
        with tempfile.TemporaryDirectory() as t:
            response=Response();response.headers['Content-Length']='999'
            with patch.object(D.urllib.request,'urlopen',return_value=response):
                self.assertEqual(D.download(Path(t),item())['status'],'failed')

    def test_redirect_to_unapproved_host_rejected(self):
        with tempfile.TemporaryDirectory() as t, patch.object(D.urllib.request,'urlopen',return_value=Response(url='https://example.com/x')):
            self.assertEqual(D.download(Path(t),item())['status'],'failed')

    def test_404_vs_ssl_and_forbidden_do_not_become_structure_absence(self):
        for error,expected in [(urllib.error.HTTPError(item()['url'],404,'Not found',{},io.BytesIO()),'unavailable'),
                               (urllib.error.HTTPError(item()['url'],403,'Forbidden',{},io.BytesIO()),'failed'),
                               (urllib.error.URLError('SSL EOF'),'failed')]:
            with tempfile.TemporaryDirectory() as t, patch.object(D.urllib.request,'urlopen',side_effect=error):
                r=D.download(Path(t),item(),attempts=1)
                self.assertEqual(r['status'],expected)
                self.assertEqual(r['url'],item()['url'])
                self.assertFalse(list(Path(t).iterdir()))

    def test_archive_is_verified_and_repeat_archive_has_no_duplicate_manifest(self):
        with tempfile.TemporaryDirectory() as t:
            project=Path(t);out=project/'results/test';out.mkdir(parents=True)
            cache=project/'data/raw/cache';(out/'summary.json').write_text('{"status":"test"}')
            with patch.object(D.urllib.request,'urlopen',return_value=Response()):r=D.download(cache,item())
            first=D.archive_files(project,out,cache,[r],[], 'a'*64)
            second=D.archive_files(project,out,cache,[r],[], 'a'*64)
            self.assertEqual(first['members'],second['members'])
            self.assertTrue((project/second['archive']).is_file())
            self.assertTrue((out/'archive_manifest.csv').is_file())

    def test_changed_binding_rejected_before_any_request(self):
        with tempfile.TemporaryDirectory() as t:
            project=Path(t);(project/'config').mkdir()
            binding=json.loads((PROJECT/'config/GSE286387_step05b_input.json').read_text())
            binding['status']='unreviewed'
            (project/'config/GSE286387_step05b_input.json').write_text(json.dumps(binding))
            with self.assertRaisesRegex(ValueError,'Unaccepted'):
                D.load_input(project)

    def test_probe_failure_leaves_bulk_pending_and_no_success(self):
        with tempfile.TemporaryDirectory() as t:
            project=Path(t)
            row=dict(protein_accession='P12345',source='PDBe_SIFTS',structure_id='1abc',chain_id='A',
                     candidate_status='candidate',species_taxid='10090',reason='pending')
            binding={'source_commit':'test'}
            with patch.object(D,'load_input',return_value=(binding,[row],[])),\
                 patch.object(D.urllib.request,'urlopen',side_effect=urllib.error.URLError('SSL')),\
                 patch.object(D.time,'sleep'):
                self.assertEqual(D.run(project,probe=True),1)
            out=project/'results/05b1_structure_download'
            self.assertFalse((out/'SUCCESS.txt').exists())
            self.assertFalse((out/'PROBE_SUCCESS.txt').exists())
            self.assertEqual(D.I.read_csv(out/'download_manifest.csv')[0]['status'],'failed')

    def test_full_acquisition_archive_and_offline_resume_in_clean_project(self):
        with tempfile.TemporaryDirectory() as t:
            project=Path(t)
            row=dict(protein_accession='P12345',source='PDBe_SIFTS',structure_id='1abc',chain_id='A',
                     candidate_status='candidate',species_taxid='10090',reason='pending')
            binding={'source_commit':'test'}
            with patch.object(D,'load_input',return_value=(binding,[row],[])), patch.object(D.urllib.request,'urlopen',return_value=Response()):
                self.assertEqual(D.run(project,archive=True),0)
            out=project/'results/05b1_structure_download'
            self.assertTrue((out/'SUCCESS.txt').exists())
            with patch.object(D,'load_input',return_value=(binding,[row],[])), patch.object(D.urllib.request,'urlopen',side_effect=AssertionError('network forbidden')):
                self.assertEqual(D.run(project,offline=True,archive=True),0)
            receipt=json.loads((out/'archive_receipt.json').read_text())
            self.assertTrue((project/receipt['archive']).is_file())


if __name__=='__main__':unittest.main()
