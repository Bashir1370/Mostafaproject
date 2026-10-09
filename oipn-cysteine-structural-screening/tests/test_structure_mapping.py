"""Strict reference mapping and SG bookkeeping precede structural quality acceptance."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

PROJECT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('mapper',PROJECT/'scripts/05b2_map_cysteines_GSE286387.py')
M=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(M)
G=M.gemmi_module()


def fixture(sequence=('ALA','CYS','GLY'), atoms=None, refs='', duplicate=False):
    poly=''.join(f'1 {n} {mon}\n' for n,mon in enumerate(sequence,1))
    if duplicate:poly+='1 2 SER\n'
    if atoms is None:atoms='1 S SG . CYS L 2 A 52 B 1 1 1 2 3 95\n'
    return G.cif.read_string('''data_test
_entity_poly.entity_id 1
_entity_poly.type 'polypeptide(L)'
_struct_asym.id L
_struct_asym.entity_id 1
loop_
_entity_poly_seq.entity_id
_entity_poly_seq.num
_entity_poly_seq.mon_id
'''+poly+refs+'''
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.pdbx_PDB_model_num
_atom_site.occupancy
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
'''+atoms).sole_block()


def candidate(source='PDBe_SIFTS'):
 return dict(protein_accession='P12345',source=source,structure_id='test',chain_id='A',sequence_start='1',sequence_end='3')


FILE=dict(file_id='id',sha256='sha',local_file='raw.cif')
REF='''
_struct_ref.id 1
_struct_ref.entity_id 1
_struct_ref.db_name UNP
_struct_ref.pdbx_db_accession P12345
_struct_ref_seq.align_id 1
_struct_ref_seq.ref_id 1
_struct_ref_seq.pdbx_strand_id A
_struct_ref_seq.seq_align_beg 1
_struct_ref_seq.seq_align_end 3
_struct_ref_seq.db_align_beg 3
_struct_ref_seq.db_align_end 5
'''


class CysMapping(unittest.TestCase):
 def test_author_numbering_insertion_and_label_chain_are_distinct(self):
  rows=M.map_candidate(fixture(),candidate(),1,'ACG',FILE)
  self.assertEqual(len(rows),1);r=rows[0]
  self.assertEqual((r['label_asym_id'],r['label_seq_id'],r['auth_seq_id'],r['insertion_code']),('L',2,'52','B'))
  self.assertEqual(r['mapping_status'],'mapped');self.assertEqual(r['structural_eligibility'],'not_assessed')

 def test_unique_subsequence_offset_and_repeated_sequence_ambiguity(self):
  r=M.map_candidate(fixture(),candidate(),1,'MMACG',FILE)[0]
  self.assertEqual(r['canonical_cys_position'],4);self.assertEqual(r['label_seq_id'],2)
  rr=M.map_candidate(fixture(),candidate(),1,'ACGACG',FILE)
  self.assertTrue(all(r['mapping_status']=='held' for r in rr))

 def test_declared_reference_segments_do_not_use_author_offsets(self):
  r=M.map_candidate(fixture(refs=REF),candidate(),1,'MMACG',FILE)[0]
  self.assertEqual((r['canonical_cys_position'],r['label_seq_id']),(4,2))
  self.assertEqual(r['mapping_method'],'checked_deposited_UNP_segments')

 def test_gapped_or_mismatched_reference_is_held_without_forced_alignment(self):
  for refs,seq in [(REF.replace('db_align_end 5','db_align_end 6'),'MMACGG'),(REF,'MMTCG')]:
   rr=M.map_candidate(fixture(refs=refs),candidate(),1,seq,FILE)
   self.assertTrue(all(r['mapping_status']=='held' for r in rr))

 def test_documented_target_substitution_is_mapped_but_held_as_changed(self):
  dif='''
_struct_ref_seq_dif.align_id 1
_struct_ref_seq_dif.seq_num 2
_struct_ref_seq_dif.pdbx_seq_db_seq_num 4
_struct_ref_seq_dif.pdbx_seq_db_accession_code P12345
_struct_ref_seq_dif.mon_id SER
_struct_ref_seq_dif.db_mon_id CYS
'''
  atoms='1 O OG . SER L 2 A 52 B 1 1 1 2 3 95\n'
  rr=M.map_candidate(fixture(sequence=('ALA','SER','GLY'),refs=REF+dif,atoms=atoms),candidate(),1,'MMACG',FILE)
  self.assertEqual(rr[0]['reason'],'TARGET_SUBSTITUTED_OR_MODIFIED')
  self.assertEqual(rr[0]['declared_changed_canonical_positions'],'4')

 def test_missing_sg_and_unobserved_residue_are_not_mapped(self):
  for atoms,reason in [('1 C CA . CYS L 2 A 52 B 1 1 1 2 3 95\n','CYS_SG_UNOBSERVED'),
                       ('1 C CA . ALA L 1 A 51 . 1 1 1 2 3 95\n','RESIDUE_UNOBSERVED')]:
   r=M.map_candidate(fixture(atoms=atoms),candidate(),1,'ACG',FILE)[0]
   self.assertEqual(r['reason'],reason);self.assertNotEqual(r['mapping_status'],'mapped')

 def test_all_models_and_occupancy_then_altloc_tie_rule(self):
  atoms='1 S SG B CYS L 2 A 52 B 1 .4 1 2 3 90\n2 S SG A CYS L 2 A 52 B 1 .6 4 5 6 95\n3 S SG . CYS L 2 A 52 B 2 1 7 8 9 98\n'
  rr=M.map_candidate(fixture(atoms=atoms),candidate(),1,'ACG',FILE)
  self.assertEqual(len(rr),2);self.assertEqual(rr[0]['selected_SG_alt_id'],'A')
  self.assertEqual(rr[1]['SG_x'],7)
  tie=atoms.replace('1 .4','1 .6');r=M.map_candidate(fixture(atoms=tie),candidate(),1,'ACG',FILE)[0]
  self.assertEqual(r['selected_SG_alt_id'],'A')

 def test_microheterogeneity_and_zero_occupancy_are_held(self):
  self.assertEqual(M.map_candidate(fixture(duplicate=True),candidate(),1,'ACG',FILE)[0]['reason'],'MICROHETEROGENEOUS_POLYMER_SEQUENCE')
  atoms='1 S SG . CYS L 2 A 52 B 1 0 1 2 3 95\n'
  self.assertEqual(M.map_candidate(fixture(atoms=atoms),candidate(),1,'ACG',FILE)[0]['reason'],'SG_ZERO_OCCUPANCY')

 def test_predicted_sequence_difference_does_not_pass(self):
  rr=M.map_candidate(fixture(),candidate('AlphaFold_DB'),1,'CCG',FILE)
  self.assertTrue(all(r['reason']=='PREDICTED_DEPOSITED_SEQUENCE_MISMATCH' for r in rr))

 def test_nopolymer_ligands_sharing_author_chain_do_not_create_candidates(self):
  atoms='1 S SG . CYS L 2 A 52 B 1 1 1 2 3 95\n2 O O . HOH W . A 99 . 1 1 4 5 6 10\n'
  self.assertEqual(len(M.map_candidate(fixture(atoms=atoms),candidate(),1,'ACG',FILE)),1)

 def test_site_audit_keeps_metadata_ambiguity_and_never_approves_quality(self):
  audit=M.site_audit({'P':'CC','Q':'C'},[],{'P'})
  self.assertEqual([r['status'] for r in audit],['held','held','unassessable'])
  self.assertTrue(all(r['structural_eligibility']=='not_assessed' for r in audit))

 def test_real_input_chain_and_all_raw_gate_reject_missing_bytes(self):
  b,seqs,candidates,links,ledger,archive,inputs=M.load_input(PROJECT)
  self.assertEqual(len(seqs),760);self.assertEqual(len(ledger),2638)
  with tempfile.TemporaryDirectory() as t:
   with self.assertRaisesRegex(ValueError,'missing or changed'):
    M.verify_raw(Path(t),ledger,archive)

 def test_representative_files_match_accepted_bytes_and_map_real_cys(self):
  b,seqs,candidates,links,ledger,archive,inputs=M.load_input(PROJECT)
  for name,expected in [('11gl',4),('AF-A0A087WRH0-F1',35)]:
   f=next(r for r in ledger if r['structure_id']==name and r['kind'].endswith('mmcif'))
   path=PROJECT/f['local_file']
   if not path.exists():continue # Unit suite remains runnable on fresh workstation; real sample validation is separately documented.
   self.assertEqual(M.file_sha(path),f['sha256']);block=G.cif.read_file(str(path)).sole_block()
   c=next(c for c in candidates if c['structure_id']==name and c['candidate_status']=='candidate')
   rr=M.map_candidate(block,c,1,seqs[c['protein_accession']],f)
   self.assertEqual(sum(r['mapping_status']=='mapped' for r in rr),expected)

 def test_end_to_end_sample_mapping_preserves_full_canonical_universe(self):
  b,seqs,candidates,links,ledger,archive,inputs=M.load_input(PROJECT)
  selected=[(n,c) for n,c in enumerate(candidates,1) if c['candidate_status']=='candidate' and c['structure_id'] in ('11gl','AF-A0A087WRH0-F1')]
  selected_candidates=[c for _,c in selected]
  _,selected_links=M.D.make_plan(selected_candidates)
  selected_ledger=[r for r in ledger if r['kind'].endswith('mmcif') and r['structure_id'] in ('11gl','AF-A0A087WRH0-F1')]
  if any(not (PROJECT/r['local_file']).exists() for r in selected_ledger):self.skipTest('Representative raw bodies not available')
  with tempfile.TemporaryDirectory() as t:
   project=Path(t)
   for r in selected_ledger:
    body=project/r['local_file'];body.parent.mkdir(parents=True,exist_ok=True);body.write_bytes((PROJECT/r['local_file']).read_bytes())
    meta=dict(file_id=r['file_id'],url=r['url'],response_url=r['response_url'],http_status=200,bytes=int(r['bytes']),sha256=r['sha256'],retrieved_utc=r['retrieved_utc'],etag=r['etag'],last_modified=r['last_modified'])
    (project/r['metadata_file']).write_text(json.dumps(meta,indent=2)+'\n')
   (project/'requirements-structure.txt').write_bytes((PROJECT/'requirements-structure.txt').read_bytes())
   with patch.object(M,'load_input',return_value=(b,seqs,selected_candidates,selected_links,selected_ledger,archive,[])):
    summary=M.run(project)
   self.assertEqual(summary['input_proteins'],760);self.assertEqual(summary['input_sites'],10799)
   self.assertEqual(summary['coordinate_files'],2);self.assertEqual(summary['sites_with_mapped_SG_options'],39)
   out=project/'results/05b2_cysteine_mapping'
   self.assertTrue((out/'SUCCESS.txt').exists())
   self.assertEqual(len(M.I.read_csv(out/'site_step_audit.csv')),10799)
   self.assertTrue(all(r['structural_eligibility']=='not_assessed' for r in M.I.read_csv(out/'site_structure_mapping.csv')))


if __name__=='__main__':unittest.main()
