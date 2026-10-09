"""Context evidence is source-bound, conformer-aware and never an automatic structural pass."""
import gzip
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import urllib.error

PROJECT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('context',PROJECT/'scripts/05b4_context_GSE286387.py')
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
OTHER = importlib.util.spec_from_file_location('quality_test',Path(__file__).with_name('test_local_quality.py'))
T = importlib.util.module_from_spec(OTHER)
OTHER.loader.exec_module(T)


def xml(pdb='11GL'):
    return gzip.compress(f'''<wwPDB-validation-information><Entry pdbid="{pdb}" PDB-revision-date="2026-07-29"/>
<ModelledSubgroup model="1" chain="A" resnum="52" icode="B" resname="CYS" altcode=" " rscc="0.9">
<clash atom="SG"/></ModelledSubgroup></wwPDB-validation-information>'''.encode())


class Response:
    status = 200
    def __init__(self,raw,url):
        self.raw,self.url = raw,url
        self.headers = {'Content-Length':str(len(raw)),'ETag':'test'}
    def __enter__(self):return self
    def __exit__(self,*args):pass
    def geturl(self):return self.url
    def read(self):return self.raw


def fixture(omit='', extra='', refs='', sequence='ACG', candidate=None):
    block, candidate_default, original, _, _ = T.model(omit=omit,sequence=sequence)
    if refs or extra:
        # Use the same explicit, complete CYS fixture but add raw partner/reference annotations.
        text=block.as_string()
        if extra:
            text += extra
        if refs:
            text += refs
        block = C.M.gemmi_module().cif.read_string(text).sole_block()
    candidate = candidate or candidate_default
    if candidate['source'] != 'AlphaFold_DB':
        original = C.M.map_candidate(block,candidate,1,sequence,T.T.FILE)
    q=C.Q.initial_row(original[0]);q['local_gate_status']='pass' if candidate['source']=='AlphaFold_DB' else 'held'
    tables={cat:list(C.M.records(block,cat)) for cat in C.ANNOTATIONS}
    atoms=C.atom_context(block,{'1'})
    return block,candidate,original[0],q,atoms,tables


class ReportAcquisition(unittest.TestCase):
    def test_official_plan_and_url_allowlist(self):
        item=C.V.plan(['11gl','11gl'])[0]
        self.assertEqual(item['url'],'https://files.wwpdb.org/pub/pdb/validation_reports/1g/11gl/11gl_validation.xml.gz')
        with self.assertRaises(ValueError):C.V.plan(['../../x'])
        with self.assertRaises(ValueError):C.V.permitted('http://files.wwpdb.org/pub/pdb/validation_reports/x')

    def test_xml_identity_and_author_insertion_matching_preserve_missing_metrics(self):
        root,entry=C.V.xml_content(xml(),'11gl');index=C.V.residue_index(root)
        status,payload=C.V.match_residue(index,'1','A','52','B','CYS','')
        self.assertEqual(status,'matched_for_manual_review')
        self.assertNotIn('rsrz',payload['attributes'])
        self.assertEqual(payload['reported_children'][0]['tag'],'clash')
        self.assertEqual(C.V.match_residue(index,'2','A','52','B','CYS','')[0],'not_reported')
        self.assertEqual(C.V.match_residue(index,'1','A','52','','CYS','')[0],'not_reported')
        with self.assertRaises(ValueError):C.V.xml_content(xml('1ABC'),'11gl')

    def test_duplicate_report_matches_are_ambiguous(self):
        root,_=C.V.xml_content(xml(),'11gl');index=C.V.residue_index(root)
        k=next(iter(index));index[k].append(index[k][0])
        self.assertEqual(C.V.match_residue(index,'1','A','52','B','CYS','')[0],'ambiguous_report_records')

    def test_cache_reuse_hashes_and_corruption(self):
        item=C.V.plan(['11gl'])[0]
        with tempfile.TemporaryDirectory() as t:
            cache=Path(t)
            with patch.object(C.V.urllib.request,'urlopen',return_value=Response(xml(),item['url'])) as request:
                first=C.V.download(cache,item,network=True)
            self.assertEqual(first['status'],'downloaded');self.assertEqual(request.call_count,1)
            with patch.object(C.V.urllib.request,'urlopen',side_effect=AssertionError('No refresh allowed')):
                second=C.V.download(cache,item,network=True)
            self.assertEqual(first['sha256'],second['sha256'])
            self.assertEqual(first['retrieved_utc'],second['retrieved_utc'])
            Path(first['local_file']).write_bytes(b'corrupt')
            self.assertEqual(C.V.download(cache,item)['status'],'failed')

    def test_offline_missing_404_and_403_are_distinct(self):
        item=C.V.plan(['11gl'])[0]
        with tempfile.TemporaryDirectory() as t:
            self.assertEqual(C.V.download(Path(t),item)['status'],'not_downloaded')
            for code,status in [(404,'unavailable'),(403,'failed')]:
                error=urllib.error.HTTPError(item['url'],code,'error',{},io.BytesIO(b'diagnostic'))
                with patch.object(C.V.urllib.request,'urlopen',side_effect=error):
                    row=C.V.download(Path(t),item,network=True)
                self.assertEqual(row['status'],status)
                self.assertEqual(row['error_excerpt'],'diagnostic')

    def test_transient_server_failure_retries_exact_url(self):
        item=C.V.plan(['11gl'])[0]
        error=urllib.error.HTTPError(item['url'],502,'bad gateway',{},io.BytesIO(b'failed'))
        with tempfile.TemporaryDirectory() as t:
            with patch.object(C.V.urllib.request,'urlopen',side_effect=[error,Response(xml(),item['url'])]) as request,patch.object(C.V.time,'sleep'):
                row=C.V.download(Path(t),item,network=True)
            self.assertEqual(row['status'],'downloaded');self.assertEqual(request.call_count,2)
            self.assertTrue(all(x.args[0].full_url==item['url'] for x in request.call_args_list))

    def test_wrong_entry_and_incomplete_transport_are_failed_not_cached(self):
        item=C.V.plan(['11gl'])[0]
        with tempfile.TemporaryDirectory() as t:
            for payload in [xml('1ABC'),b'not-gzip']:
                with patch.object(C.V.urllib.request,'urlopen',return_value=Response(payload,item['url'])):
                    row=C.V.download(Path(t),item,network=True)
                self.assertEqual(row['status'],'failed')
            response=Response(xml(),item['url']);response.headers['Content-Length']=str(len(response.raw)+1)
            with patch.object(C.V.urllib.request,'urlopen',return_value=response):
                self.assertEqual(C.V.download(Path(t),item,network=True)['status'],'failed')
            self.assertFalse(list(Path(t).glob('*.metadata.json')))


class ContextEvidence(unittest.TestCase):
    def test_target_and_neighbor_gaps_are_evidence_not_pass(self):
        b,c,m,q,a,t=fixture(omit='CB')
        row,_=C.evaluate_option(b,c,'ACG',m,q,a,t)
        self.assertEqual(row['missing_coherent_target_atoms'],'CB')
        self.assertTrue(json.loads(row['neighbor_atom_gaps']))
        self.assertEqual(row['structural_eligibility'],'held')

    def test_alternate_target_atoms_are_not_union_completed(self):
        b,c,m,q,_,t=fixture()
        table=b.find_mmcif_category('_atom_site.')
        names=list(table.tags)
        for r in table:
            if r[names.index('_atom_site.label_atom_id')]=='CB':r[names.index('_atom_site.label_alt_id')]='A'
        a=C.atom_context(b,{'1'})
        row,_=C.evaluate_option(b,c,'ACG',m,q,a,t)
        self.assertEqual(row['missing_coherent_target_atoms'],'CB')

    def test_assembly_recipes_preserved_without_generation_or_selection(self):
        extra='''\n_pdbx_struct_assembly.id 1
_pdbx_struct_assembly.details author_defined
_pdbx_struct_assembly_gen.assembly_id 1
_pdbx_struct_assembly_gen.oper_expression '(1-3)(4,5)'
_pdbx_struct_assembly_gen.asym_id_list L
'''
        b,c,m,q,a,t=fixture(extra=extra)
        row,_=C.evaluate_option(b,c,'ACG',m,q,a,t)
        self.assertEqual(row['candidate_assembly_ids'],'1')
        self.assertEqual(row['assembly_definition_status'],'deposited_recipes_not_generated')
        self.assertEqual(t['_pdbx_struct_assembly_gen.'][0]['oper_expression'],'(1-3)(4,5)')

    def test_other_subchains_ligands_and_unknown_coordinate_context_are_retained(self):
        b,c,m,q,_,t=fixture()
        # Add a separate nonpolymer atom without altering canonical chain correspondence.
        table=b.find_mmcif_category('_atom_site.')
        data={str(k).split('.',1)[1]:'.' for k in table.tags}
        data.update(id='99',type_symbol='ZN',label_atom_id='ZN',label_comp_id='ZN',label_asym_id='Z',
                    auth_asym_id='A',auth_seq_id='500',pdbx_PDB_model_num='1',occupancy='1',Cartn_x='1',Cartn_y='0',Cartn_z='0')
        table.append_row([data[str(k).split('.',1)[1]] for k in table.tags])
        a=C.atom_context(b,{'1'})
        row,_=C.evaluate_option(b,c,'ACG',m,q,a,t)
        self.assertEqual(row['neighbor_other_subchains'],'Z')
        self.assertEqual(row['neighbor_nonpolymer_components'],'ZN')
        self.assertTrue(json.loads(row['neighbor_unrecognized_components']))
        table[-1][list(table.tags).index('_atom_site.Cartn_x')]='?'
        a=C.atom_context(b,{'1'})
        row,_=C.evaluate_option(b,c,'ACG',m,q,a,t)
        self.assertEqual(row['reason_code'],'MODEL_COORDINATE_CONTEXT_UNRESOLVED')

    def test_failed_local_option_remains_unassessed_without_reintroducing_it(self):
        b,c,m,q,_,t=fixture();q['local_gate_status']='excluded'
        row,_=C.evaluate_option(b,c,'ACG',m,q,None,t)
        self.assertEqual(row['context_status'],'not_assessed')
        self.assertEqual(row['local_gate_status'],'excluded')

    def test_declared_neighbor_mutation_uses_canonical_mapping_not_author_number(self):
        dif='''
_struct_ref_seq_dif.align_id 1
_struct_ref_seq_dif.seq_num 1
_struct_ref_seq_dif.pdbx_seq_db_seq_num 3
_struct_ref_seq_dif.pdbx_seq_db_accession_code P12345
_struct_ref_seq_dif.mon_id ALA
_struct_ref_seq_dif.db_mon_id SER
'''
        b,c,m,q,a,t=fixture(refs=T.T.REF+dif,sequence='MMSCG',candidate=T.T.candidate())
        row,_=C.evaluate_option(b,c,'MMSCG',m,q,a,t)
        self.assertEqual(row['declared_mutations_in_neighborhood'],'3')
        self.assertEqual(m['canonical_cys_position'],4)
        self.assertEqual(row['structural_eligibility'],'held')

    def test_sg_association_mismatch_fails(self):
        b,c,m,q,a,t=fixture();m['SG_x']=1
        with self.assertRaisesRegex(ValueError,'SG/context identity mismatch'):
            C.evaluate_option(b,c,'ACG',m,q,a,t)

    def test_frozen_actual_input_and_full_original_provenance(self):
        values=C.load_input(PROJECT)
        self.assertEqual(len(values[5]),14310)
        self.assertEqual(len(values[6]),10799)
        self.assertEqual(len(C.V.plan(r['structure_id'] for r in values[3] if r['kind']=='pdb_mmcif')),445)

    def test_real_reference_sample_context_and_brute_force_neighbor_geometry(self):
        _,seqs,cs,ledger,_,quality,_,originals,_=C.load_input(PROJECT)
        raw=next(r for r in ledger if r['structure_id']=='11gl' and r['kind']=='pdb_mmcif')
        if not (PROJECT/raw['local_file']).is_file():self.skipTest('Representative raw coordinates unavailable')
        b=C.M.gemmi_module().cif.read_file(str(PROJECT/raw['local_file'])).sole_block()
        self.assertEqual(C.M.file_sha(PROJECT/raw['local_file']),raw['sha256'])
        tables={category:list(C.M.records(b,category)) for category in C.ANNOTATIONS}
        atoms=C.atom_context(b,{'1'})
        options=[r for r in quality if r['structure_id']=='11gl']
        self.assertEqual(len(options),8)
        for q in options:
            m=originals[(q['site_id'],q['candidate_row'],q['label_asym_id'],q['model_id'])]
            row,_=C.evaluate_option(b,cs[int(q['candidate_row'])-1],seqs[q['protein_accession']],m,q,atoms,tables)
            target=atoms[2][('1',m['selected_SG_atom_id'])]
            near={C.atom_key(a) for a in atoms[2].values() if a['type_symbol'].upper() not in ('H','D') and
                  C.atom_key(a)!=C.atom_key(target) and a['occupancy']!='0' and
                  sum((x-y)**2 for x,y in zip(a['xyz'],target['xyz']))<=36}
            self.assertEqual(row['neighbor_residue_count'],len(near))
            self.assertEqual(row['neighborhood_frame'],'deposited_ASU')
            self.assertEqual(row['structural_eligibility'],'held')

    def test_end_to_end_context_keeps_full_site_audit_and_invalidates_stale_success(self):
        b,seqs,cs,ledger,archive,quality,sites,originals,_=C.load_input(PROJECT)
        selected=[r for r in quality if r['structure_id'] in ('11gl','AF-A0A087WRH0-F1')]
        small=[r for r in ledger if r['structure_id'] in ('11gl','AF-A0A087WRH0-F1')]
        if any(not (PROJECT/r['local_file']).is_file() for r in small):self.skipTest('Representative raw coordinates unavailable')
        with tempfile.TemporaryDirectory() as t:
            project=Path(t)
            for r in small:
                f=project/r['local_file'];f.parent.mkdir(parents=True,exist_ok=True);f.write_bytes((PROJECT/r['local_file']).read_bytes())
                meta=dict(file_id=r['file_id'],url=r['url'],response_url=r['response_url'],http_status=200,bytes=int(r['bytes']),
                          sha256=r['sha256'],retrieved_utc=r['retrieved_utc'],etag=r['etag'],last_modified=r['last_modified'])
                (project/r['metadata_file']).write_text(json.dumps(meta,indent=2)+'\n')
            # When present, replay the independently acquired official report through the real cache gate.
            probe=PROJECT/'results/05b4_structure_context/probe/validation_download_manifest.csv'
            report_cached=False
            if probe.is_file():
                records=C.I.read_csv(probe)
                report=next((r for r in records if r['structure_id']=='11gl' and r['status']=='downloaded'),None)
                if report and (PROJECT/report['local_file']).is_file():
                    planned=C.V.plan(['11gl']);plan_hash=C.V.sha(json.dumps(planned,sort_keys=True).encode())
                    cache=project/'data/raw/experimental_validation/step05b4'/plan_hash;cache.mkdir(parents=True)
                    for field in ('local_file','metadata_file'):
                        (cache/Path(report[field]).name).write_bytes((PROJECT/report[field]).read_bytes())
                    report_cached=True
            values=b,seqs,cs,small,archive,selected,sites,originals,[]
            with patch.object(C,'load_input',return_value=values):
                summary=C.run(project)
            out=project/'results/05b4_structure_context'
            self.assertEqual(summary['mapped_options'],43)
            self.assertEqual(summary['input_sites'],10799)
            self.assertEqual(summary['coordinate_files'],2)
            self.assertEqual(summary['final_structural_passes'],0)
            if report_cached:
                self.assertEqual(summary['experimental_target_validation_status_counts'],{'matched_for_manual_review':8})
            self.assertEqual(summary['inherited_site_status_counts'],{'excluded':3042,'held':6964,'unassessable':793})
            self.assertTrue((out/'SUCCESS.txt').is_file())
            self.assertTrue(all(r['review_status']=='pending' for r in C.I.read_csv(out/'manual_review_queue.csv')))
            (project/small[0]['local_file']).write_bytes(b'bad')
            with patch.object(C,'load_input',return_value=values):
                with self.assertRaises(ValueError):C.run(project)
            self.assertFalse((out/'SUCCESS.txt').exists())


if __name__ == '__main__':unittest.main()
