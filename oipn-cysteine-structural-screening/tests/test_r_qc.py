"""Integration tests with actual audited data, when R/DESeq2/jsonlite are available."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / 'scripts/01b_qc_GSE286387.R'
INPUT = PROJECT / 'results/01_dataset_audit'
HAS_INPUT = (INPUT / 'count_matrix.tsv').is_file()
HAS_R = shutil.which('Rscript') is not None


@unittest.skipUnless(HAS_INPUT and HAS_R, 'Requires actual Step 01 inputs and Rscript')
class SampleQC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        r = subprocess.run(['Rscript', '-e',
                            "quit(status=if(all(vapply(c('DESeq2','jsonlite'),requireNamespace,logical(1),quietly=TRUE))) 0 else 1)"],
                           capture_output=True, timeout=120)
        if r.returncode:
            raise unittest.SkipTest('Requires installed DESeq2 and jsonlite')

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        for folder in ['config', 'docs/audits', 'results/01_dataset_audit']:
            shutil.copytree(PROJECT / folder, self.root / folder)
        (self.root / 'data').mkdir()
        (self.root / 'data/raw').symlink_to(PROJECT / 'data/raw', target_is_directory=True)

    def tearDown(self):
        self.tmp.cleanup()

    def run_script(self, via_source=False):
        env = dict(os.environ, OIPN_PROJECT_DIR=str(self.root))
        args = ['Rscript', str(SCRIPT)] if not via_source else ['Rscript', '-e', 'source(' + repr(str(SCRIPT)) + ')']
        return subprocess.run(args, env=env, capture_output=True, text=True, timeout=240)

    def test_real_inputs_create_qc_without_approving_samples(self):
        run = self.run_script(via_source=True)
        self.assertEqual(run.returncode, 0, run.stderr)
        out = self.root / 'results/01b_sample_qc'
        self.assertTrue((out / 'SUCCESS.txt').is_file())
        self.assertFalse((out / 'FAILURE.txt').exists())
        for name in ['pca_top500.png', 'sample_distance.png', 'sample_correlation.png',
                     'library_size.png', 'QC_plots.pdf', 'sessionInfo.txt', 'input_checksums.csv']:
            self.assertGreater((out / name).stat().st_size, 100)
        with (out / 'sample_metrics.csv').open() as f:
            rows = list(csv.DictReader(f))
        self.assertEqual(len(rows), 10)
        self.assertEqual({r['inclusion_status'] for r in rows}, {'held'})
        self.assertIn('18538', (out / 'qc_report.md').read_text())

    def test_changed_count_is_rejected(self):
        path = self.root / 'results/01_dataset_audit/count_matrix.tsv'
        lines = path.read_text().splitlines()
        row = lines[1].split('\t')
        row[1] = str(int(row[1]) + 1)
        lines[1] = '\t'.join(row)
        path.write_text('\n'.join(lines) + '\n')
        run = self.run_script()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('differs from raw input', run.stderr)
        self.assertFalse((self.root / 'results/01b_sample_qc/SUCCESS.txt').exists())

    def test_changed_condition_is_rejected(self):
        path = self.root / 'results/01_dataset_audit/sample_manifest.csv'
        with path.open() as f:
            reader = csv.DictReader(f)
            fields, rows = reader.fieldnames, list(reader)
        rows[0]['condition'] = 'oxaliplatin' if rows[0]['condition'] == 'control' else 'control'
        with path.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        run = self.run_script()
        self.assertNotEqual(run.returncode, 0)
        self.assertIn('differs from audited mapping', run.stderr)


if __name__ == '__main__':
    unittest.main()
