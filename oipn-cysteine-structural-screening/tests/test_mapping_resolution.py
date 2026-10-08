"""Evidence-driven representative selection must not become a reviewed-entry shortcut."""
import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

PROJECT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('resolver', PROJECT / 'scripts/03b_resolve_GSE286387.py')
R = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(R)
G = 'ENSMUSG00000000001'


def entry(acc):
    return {'primaryAccession': acc, 'organism': {'taxonId': 10090},
            'sequence': {'value': 'MACK', 'length': 4},
            'uniProtKBCrossReferences': [{'database': 'Ensembl', 'id': 'ENSMUST00000000001.2',
                                        'properties': [{'key': 'GeneId', 'value': G + '.4'}]}]}


def group(parent='P12345', members=('P12345', 'A0A1234567'), gene=G):
    return {'gene': {'accession': parent, 'geneNameType': 'MOD', 'geneName': 'MGI:0'},
            'relatedGene': [{'accession': acc, 'geneNameType': 'Ensembl', 'geneName': gene} for acc in members]}


class RepresentativeEvidence(unittest.TestCase):
    def setUp(self):
        self.options = [entry('P12345'), entry('A0A1234567')]

    def test_explicit_group_resolves_and_duplicate_evidence_is_harmless(self):
        for groups in [[group()], [group(), group()]]:
            choice, reason, evidence = R.choose_representative(G, self.options, groups)
            self.assertEqual(choice[0]['primaryAccession'], 'P12345')
            self.assertEqual(reason, 'RESOLVED_EXPLICIT_GENE_CENTRIC_REPRESENTATIVE')

    def test_review_status_alone_and_gene_symbol_never_resolve(self):
        self.options[0]['entryType'] = 'UniProtKB reviewed (Swiss-Prot)'
        self.options[1]['entryType'] = 'UniProtKB unreviewed (TrEMBL)'
        self.assertIsNone(R.choose_representative(G, self.options, [])[0])
        self.assertIsNone(R.choose_representative(G, self.options, [group(gene='MyGeneSymbol')])[0])

    def test_unaccounted_candidate_and_competing_parents_remain_held(self):
        self.assertEqual(R.choose_representative(G, self.options, [group(members=('P12345',))])[1],
                         'CANDIDATES_NOT_FULLY_ACCOUNTED_FOR')
        self.assertEqual(R.choose_representative(G, self.options, [group(), group(parent='A0A1234567')])[1],
                         'MULTIPLE_GENE_CENTRIC_REPRESENTATIVES')

    def test_parent_without_original_gene_link_or_wrong_isoform_remains_held(self):
        self.assertEqual(R.choose_representative(G, self.options, [group(parent='Q12345')])[1],
                         'REPRESENTATIVE_LACKS_ORIGINAL_EXACT_GENE_LINK')
        self.assertEqual(R.choose_representative(G, self.options, [group(parent='P12345-2')])[1],
                         'REPRESENTATIVE_IS_NOT_CURRENT_DISPLAYED_ISOFORM')

    def test_representative_sequence_still_requires_full_quality(self):
        self.options[0]['proteinDescription'] = {'flag': 'Fragment'}
        self.assertEqual(R.choose_representative(G, self.options, [group()])[1], 'REPRESENTATIVE_FRAGMENT_SEQUENCE')

    def test_modified_baseline_input_is_rejected_before_resolution(self):
        with tempfile.TemporaryDirectory() as temp:
            project = Path(temp); baseline = project / 'results/03_protein_mapping'
            baseline.mkdir(parents=True)
            (baseline / 'SUCCESS.txt').write_text('STEP03_GENERATED_REVIEW_PENDING\n')
            (project / 'input.txt').write_text('modified')
            R.M.write_csv(baseline / 'input_checksums.csv',
                          [{'input': 'input.txt', 'sha256': R.M.sha(b'original')}], ['input', 'sha256'])
            with self.assertRaisesRegex(ValueError, 'checksum changed'):
                R.load_baseline(project, {}, [])


@unittest.skipUnless((PROJECT / 'results/03b_mapping_resolution/SUCCESS.txt').exists(), 'Run Step 03b first')
class ResolvedOutput(unittest.TestCase):
    def test_baseline_passes_preserved_and_all_final_fasta_hashes_match(self):
        base = R.M.read_csv(PROJECT / 'results/03_protein_mapping/gene_protein_mapping.csv')
        out = PROJECT / 'results/03b_mapping_resolution'
        final = R.M.read_csv(out / 'gene_protein_mapping.csv')
        indexed = {r['stable_gene_id']: r for r in final}
        self.assertEqual(len(final), 854)
        self.assertEqual({r['stable_gene_id'] for r in base}, set(indexed))
        for row in base:
            if row['mapping_status'] == 'pass':
                self.assertEqual(indexed[row['stable_gene_id']]['protein_accession'], row['protein_accession'])
                self.assertEqual(indexed[row['stable_gene_id']]['sequence_checksum'], row['sequence_checksum'])
        seqs, acc = {}, None
        for line in (out / 'canonical_sequences.fasta').read_text().splitlines():
            if line.startswith('>'):
                acc = line[1:].split()[0]; self.assertNotIn(acc, seqs); seqs[acc] = ''
            else:
                seqs[acc] += line
        passed = [r for r in final if r['mapping_status'] == 'pass']
        self.assertEqual(set(seqs), {r['protein_accession'] for r in passed})
        for row in passed:
            self.assertEqual(R.M.sha(seqs[row['protein_accession']].encode()), row['sequence_checksum'])
        evidence = R.M.read_csv(out / 'gene_centric_evidence.csv')
        for row in passed:
            if row['resolution_method'] == 'explicit_gene_centric_evidence':
                self.assertTrue(any(e['stable_gene_id'] == row['stable_gene_id'] and
                                    e['covers_all_original_candidates'] == 'True' and
                                    R.base_accession(e['representative_identifier']) == row['protein_accession']
                                    for e in evidence))


if __name__ == '__main__':
    unittest.main()
