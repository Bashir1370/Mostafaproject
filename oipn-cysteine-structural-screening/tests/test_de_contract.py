"""Independent checks of the real Step 02 result and frozen-input gates."""
import csv
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / 'results/02_differential_expression'


def rows(name):
    with (OUTPUT / name).open() as handle:
        return list(csv.DictReader(handle))


@unittest.skipUnless((OUTPUT / 'SUCCESS.txt').exists(), 'Run Step 02 on audited inputs first')
class DiscoveryContract(unittest.TestCase):
    def test_bh_uses_all_valid_biotypes(self):
        valid = [r for r in rows('all_genes_de.csv') if r['valid_p_for_BH'] == 'TRUE']
        self.assertGreater(len({r['biotype'] for r in valid}), 1)
        ordered = sorted(valid, key=lambda r: float(r['pvalue']))
        adjusted = 1.0
        for rank in range(len(ordered), 0, -1):
            r = ordered[rank - 1]
            adjusted = min(adjusted, float(r['pvalue']) * len(ordered) / rank)
            self.assertTrue(math.isclose(float(r['padj']), adjusted, rel_tol=1e-10, abs_tol=1e-14))

    def test_complete_audit_and_unrestricted_significant_coding_selection(self):
        all_rows = rows('all_genes_de.csv')
        ids = {r['stable_gene_id'] for r in all_rows}
        self.assertEqual(len(ids), len(all_rows))
        self.assertEqual(ids, {r['object_id'] for r in rows('step_audit.csv')})
        self.assertEqual({r['stable_gene_id'] for r in rows('tested_gene_universe.csv')},
                         {r['stable_gene_id'] for r in all_rows if r['prefilter_status'] == 'pass'})
        expected = {r['stable_gene_id'] for r in all_rows
                    if r['padj'] != 'NA' and float(r['padj']) < .05 and r['biotype'] == 'protein_coding'}
        discovery = rows('significant_protein_coding_degs.csv')
        self.assertEqual(expected, {r['stable_gene_id'] for r in discovery})
        self.assertEqual({'up', 'down'}, {r['direction'] for r in discovery})
        self.assertTrue(any(abs(float(r['log2FC'])) < 1 for r in discovery))


@unittest.skipUnless(shutil.which('Rscript') and (OUTPUT / 'SUCCESS.txt').exists(),
                     'Requires Rscript and completed Step 02 inputs')
class FrozenInputGates(unittest.TestCase):
    def check_rejected(self, relative_path, message):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            for folder in ['config', 'docs/audits', 'results/01_dataset_audit', 'results/01b_sample_qc']:
                shutil.copytree(PROJECT / folder, root / folder)
            path = root / relative_path
            # Even a formatting-only change must invalidate frozen provenance.
            with path.open('a') as handle:
                handle.write('\n')
            run = subprocess.run(['Rscript', str(PROJECT / 'scripts/02_DE_GSE286387.R')],
                                 env=dict(os.environ, OIPN_PROJECT_DIR=str(root)),
                                 capture_output=True, text=True, timeout=120)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn(message, run.stderr)
            self.assertFalse((root / 'results/02_differential_expression/SUCCESS.txt').exists())
            self.assertTrue((root / 'results/02_differential_expression/FAILURE.txt').exists())

    def test_modified_count_rejected_before_fit(self):
        self.check_rejected('results/01_dataset_audit/count_matrix.tsv', 'Count matrix changed')

    def test_modified_sample_manifest_rejected_before_fit(self):
        self.check_rejected('config/GSE286387_samples.csv', 'Frozen sample manifest was modified')


if __name__ == '__main__':
    unittest.main()
