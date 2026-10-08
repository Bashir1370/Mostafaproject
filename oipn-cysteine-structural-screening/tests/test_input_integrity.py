"""Regression checks for count inflation, fractional input and sample misassignment."""
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('audit', Path(__file__).resolve().parents[1] / 'scripts/01_audit_GSE286387.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class InputIntegrity(unittest.TestCase):
    def read(self, rows, titles=('sample1', 'sample2')):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d) / 'counts.gz'
            with gzip.open(path, 'wt') as f:
                f.write('\t'.join(audit.ANNOTATIONS + ['sample1', 'sample2']) + '\n')
                f.write(rows)
            return audit.read_counts(path, titles)

    def test_annotation_duplicates_do_not_double_counts(self):
        _, genes, annotations, duplicates = self.read(
            'ENSMUSG00000000001\t1\tGene\tprotein_coding\tMGI\t10\t20\n'
            'ENSMUSG00000000001\t2\tGene\tprotein_coding\tMGI\t10\t20\n')
        self.assertEqual(genes['ENSMUSG00000000001']['values'], (10, 20))
        self.assertEqual(len(annotations), 2)
        self.assertEqual(len(duplicates), 1)
        self.assertEqual(genes['ENSMUSG00000000001']['entrez'], {'1', '2'})

    def test_conflicting_duplicate_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'Conflicting duplicate'):
            self.read('ENSMUSG1\t1\tGene\tprotein_coding\tMGI\t10\t20\n'
                      'ENSMUSG1\t2\tGene\tprotein_coding\tMGI\t11\t20\n')

    def test_fractional_negative_missing_nonfinite_rejected(self):
        for value in ['1.5', '-1', 'NaN', 'Infinity', '']:
            with self.subTest(value=value), self.assertRaises(Exception):
                audit.integer(value)

    def test_sample_mismatch_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'exactly match GEO'):
            self.read('ENSMUSG1\t1\tGene\tprotein_coding\tMGI\t10\t20\n', ('wrong', 'sample2'))


if __name__ == '__main__':
    unittest.main()
