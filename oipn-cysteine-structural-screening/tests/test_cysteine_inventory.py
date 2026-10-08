"""Check residue numbering, provenance gates and unique protein/site accounting."""
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

PROJECT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('inventory', PROJECT / 'scripts/04_inventory_GSE286387.py')
I = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(I)
G1, G2 = 'ENSMUSG00000000001', 'ENSMUSG00000000002'


class CysteineInventory(unittest.TestCase):
    def test_numbering_termini_adjacent_cys_and_selenocysteine(self):
        proteins, sites = I.inventory({'P1': 'CCUAC', 'P2': 'MUA'},
                                      {'P1': {G1}, 'P2': {G2}}, 'test', '0.1.0')
        self.assertEqual([s['canonical_cys_position'] for s in sites], [1, 2, 5])
        self.assertEqual([s['site_id'] for s in sites], ['P1:C1', 'P1:C2', 'P1:C5'])
        self.assertEqual(proteins[0]['cys_percent'], 60.0)
        self.assertEqual(proteins[1]['inventory_status'], 'excluded')
        self.assertEqual(proteins[1]['n_cys_total'], 0)

    def test_multiple_genes_do_not_duplicate_protein_or_sites(self):
        proteins, sites = I.inventory({'P1': 'MCC'}, {'P1': {G2, G1}}, 'test', '0.1.0')
        self.assertEqual(len(proteins), 1)
        self.assertEqual(len(sites), 2)
        self.assertEqual(proteins[0]['stable_gene_ids'], G1 + ';' + G2)

    def test_fasta_rejects_ambiguous_duplicate_and_wrong_species(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / 'input.fasta'
            for text in [f'>P1 taxid=10090 genes={G1}\nMCX\n',
                         f'>P1 taxid=10090 genes={G1}\nMC\n>P1 taxid=10090 genes={G1}\nMC\n',
                         f'>P1 taxid=9606 genes={G1}\nMC\n',
                         f'>P1 taxid=10090 genes={G1}\n']:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    I.read_fasta(path)

    def test_changed_frozen_input_removes_success_and_records_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            (project / 'config').mkdir()
            binding = json.loads((PROJECT / 'config/GSE286387_step04_input.json').read_text())
            for relative in binding['input_sha256']:
                target = project / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(PROJECT / relative, target)
            (project / 'config/GSE286387_step04_input.json').write_text(json.dumps(binding))
            fasta = project / binding['fasta_file']
            fasta.write_bytes(fasta.read_bytes() + b'\n')
            out = project / 'results/04_cysteine_inventory'
            out.mkdir(parents=True)
            (out / 'SUCCESS.txt').write_text('old success')
            with self.assertRaisesRegex(ValueError, 'Frozen input checksum changed'):
                I.run(project)
            self.assertFalse((out / 'SUCCESS.txt').exists())
            self.assertTrue((out / 'FAILURE.txt').exists())
            self.assertFalse(list(out.glob('inventory_stage_*')))

    def test_real_snapshot_inventory_and_gene_coverage(self):
        binding, rows, sequences, associated, inputs = I.load_input(PROJECT)
        proteins, sites = I.inventory(sequences, associated, binding['uniprot_sequence_release'], '0.1.0')
        self.assertEqual(len(rows), 854)
        self.assertEqual(len(proteins), 787)
        expected = {(a, m.start() + 1) for a, seq in sequences.items() for m in re.finditer('C', seq)}
        actual = {(r['protein_accession'], r['canonical_cys_position']) for r in sites}
        self.assertEqual(actual, expected)
        self.assertEqual(len(sites), len(actual))
        self.assertEqual(len(sites), 10799)
        self.assertEqual(sum(r['n_cys_total'] == 0 for r in proteins), 27)
        self.assertEqual(sum(r['n_cys_total'] > 0 for r in proteins), 760)
        # A clean checkout needs no prior results directory or raw network cache.
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp)
            paths = list(binding['input_sha256']) + ['config/GSE286387_step04_input.json',
                                                     'scripts/04_inventory_GSE286387.py']
            for relative in paths:
                target = project / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(PROJECT / relative, target)
            result = subprocess.run([sys.executable, str(project / 'scripts/04_inventory_GSE286387.py')],
                                    capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            output = project / 'results/04_cysteine_inventory'
            self.assertTrue((output / 'SUCCESS.txt').exists())
            self.assertFalse((output / 'FAILURE.txt').exists())
            gene_audit = I.read_csv(output / 'gene_step_audit.csv')
            self.assertEqual(len(gene_audit), 854)
            self.assertEqual(len({r['object_id'] for r in gene_audit}), 854)
            self.assertEqual(sum(r['status'] == 'held' for r in gene_audit), 64)
            self.assertEqual(sum(r['status'] == 'unassessable' for r in gene_audit), 1)
            self.assertEqual(len(I.read_csv(output / 'cysteine_inventory.csv')), 10799)


if __name__ == '__main__':
    unittest.main()
