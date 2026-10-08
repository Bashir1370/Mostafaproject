"""Scientific identity/ambiguity and snapshot-integrity tests for Step 03."""
import csv
import copy
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest
from unittest.mock import patch

PROJECT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('mapping', PROJECT / 'scripts/03_map_GSE286387.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
G = 'ENSMUSG00000000001'


def entry(acc='P12345', taxid=10090, isoform=None):
    ref = {'database': 'Ensembl', 'id': 'ENSMUST00000000001.4',
           'properties': [{'key': 'GeneId', 'value': G + '.7'},
                          {'key': 'ProteinId', 'value': 'ENSMUSP00000000001.3'}]}
    if isoform:
        ref['isoformId'] = isoform
    return {'primaryAccession': acc, 'organism': {'taxonId': taxid},
            'entryType': 'UniProtKB reviewed (Swiss-Prot)',
            'sequence': {'value': 'MACK', 'length': 4}, 'uniProtKBCrossReferences': [ref]}


class MappingIdentity(unittest.TestCase):
    def test_versioned_gene_and_explicit_canonical_isoform(self):
        e = entry(isoform='P12345-2')
        e['comments'] = [{'commentType': 'ALTERNATIVE PRODUCTS', 'isoforms': [
            {'isoformSequenceStatus': 'Displayed', 'isoformIds': ['P12345-2']}]}]
        self.assertIn(G, M.gene_refs(e))
        result = M.classify(G, [e])
        self.assertEqual(result[0], 'pass')
        self.assertEqual(result[3], 'P12345-2')  # Canonical is not assumed to be isoform 1.

    def test_reviewed_and_unreviewed_candidates_remain_ambiguous(self):
        other = entry('A0A1234567')
        other['entryType'] = 'UniProtKB unreviewed (TrEMBL)'
        for options in [[entry(), other], [other, entry()]]:
            self.assertEqual(M.classify(G, options)[:2], ('held', 'MULTIPLE_MOUSE_UNIPROT_ENTRIES'))

    def test_noncanonical_only_and_wrong_species_are_not_accepted(self):
        e = entry(isoform='P12345-2')
        e['comments'] = [{'commentType': 'ALTERNATIVE PRODUCTS', 'isoforms': [
            {'isoformSequenceStatus': 'Displayed', 'isoformIds': ['P12345-1']}]}]
        self.assertEqual(M.classify(G, [e])[1], 'ONLY_NONCANONICAL_OR_UNRESOLVED_ISOFORM_XREF')
        self.assertEqual(M.classify(G, [entry(taxid=9606)])[1], 'WRONG_SPECIES')
        self.assertEqual(M.classify(G, [])[1], 'NO_CURRENT_EXACT_ENSEMBL_XREF')

    def test_sequence_ambiguity_fragment_and_length_are_unassessable(self):
        for mutation, reason in [({'sequence': {'value': 'MXCK', 'length': 4}}, 'AMBIGUOUS_OR_INVALID_SEQUENCE'),
                                 ({'sequence': {'value': 'MACK', 'length': 3}}, 'SEQUENCE_LENGTH_MISMATCH'),
                                 ({'proteinDescription': {'flag': 'Fragment'}}, 'FRAGMENT_SEQUENCE')]:
            e = entry(); e.update(mutation)
            self.assertEqual(M.classify(G, [e])[:2], ('unassessable', reason))
        e = entry(); e['sequence'] = {'value': 'MAUO', 'length': 4}
        self.assertEqual(M.classify(G, [e])[0], 'pass')  # U/O are valid; Cys presence is Step 04.


class ReferenceIntegrity(unittest.TestCase):
    def test_corrupt_cached_response_is_rejected(self):
        raw = b'{"results":[]}'
        meta = {'url': 'url', 'sha256': M.sha(raw), 'release': 'test', 'release_date': 'test'}
        M.validated_page(raw, meta, 'url')
        with self.assertRaisesRegex(ValueError, 'Corrupt'):
            M.validated_page(raw + b' ', meta, 'url')

    def test_pagination_completion_and_missing_page(self):
        pages = [({'results': [entry()]}, {'total_results': 2, 'next_url': 'next'}),
                 ({'results': [entry('Q12345')]}, {'total_results': 2, 'next_url': ''})]
        with patch.object(M, 'get_page', side_effect=copy.deepcopy(pages)):
            result, provenance = M.fetch_batch([G], Path('.'), 1)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(provenance), 2)
        pages[0][1]['next_url'] = ''
        with patch.object(M, 'get_page', side_effect=pages):
            with self.assertRaisesRegex(ValueError, 'Incomplete UniProt pagination'):
                M.fetch_batch([G], Path('.'), 1)

    @unittest.skipUnless((PROJECT / 'config/GSE286387_step03_input.json').exists(), 'Requires accepted binding')
    def test_modified_discovery_is_rejected_before_network(self):
        with tempfile.TemporaryDirectory() as temp:
            p = Path(temp)
            shutil.copytree(PROJECT / 'config', p / 'config')
            binding = json.loads((p / 'config/GSE286387_step03_input.json').read_text())
            src = PROJECT / binding['discovery_file']
            if not src.exists():
                self.skipTest('Requires committed workstation discovery file')
            dst = p / binding['discovery_file']; dst.parent.mkdir(parents=True)
            dst.write_bytes(src.read_bytes() + b'\n')
            with self.assertRaisesRegex(ValueError, 'discovery CSV changed'):
                M.load_input(p)


@unittest.skipUnless((PROJECT / 'results/03_protein_mapping/SUCCESS.txt').exists(), 'Run Step 03 first')
class RealOutputContract(unittest.TestCase):
    def test_all_genes_audited_and_fasta_matches_accepted_mapping(self):
        out = PROJECT / 'results/03_protein_mapping'
        rows = M.read_csv(out / 'gene_protein_mapping.csv')
        audit = M.read_csv(out / 'step_audit.csv')
        self.assertEqual(len(rows), 854)
        self.assertEqual({r['stable_gene_id'] for r in rows}, {r['object_id'] for r in audit})
        seqs, accession = {}, None
        for line in (out / 'canonical_sequences.fasta').read_text().splitlines():
            if line.startswith('>'):
                accession = line[1:].split()[0]; self.assertNotIn(accession, seqs); seqs[accession] = ''
            else:
                seqs[accession] += line
        passed = [r for r in rows if r['mapping_status'] == 'pass']
        self.assertEqual(set(seqs), {r['protein_accession'] for r in passed})
        for row in passed:
            seq = seqs[row['protein_accession']]
            self.assertEqual(M.sha(seq.encode()), row['sequence_checksum'])
            self.assertEqual(len(seq), int(row['sequence_length']))
            self.assertEqual(row['species_taxid'], '10090')
        for row in rows:
            if row['mapping_status'] != 'pass':
                self.assertEqual(row['protein_accession'], '')


if __name__ == '__main__':
    unittest.main()
