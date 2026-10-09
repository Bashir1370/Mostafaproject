#!/usr/bin/env python3
"""Offline, hash-gated canonical Cys mapping; local quality/assembly approval is later work."""
import argparse
from collections import defaultdict, Counter
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import shutil
import sys
import tempfile

SPEC = importlib.util.spec_from_file_location('download', Path(__file__).with_name('05b1_download_GSE286387.py'))
D = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(D)
I = D.I
VERSION = '1.0.0'
GEMMI_VERSION = '0.7.5'


class MappingHold(ValueError):
    pass


def need(condition, reason):
    if not condition:
        raise MappingHold(reason)


def gemmi_module():
    try:
        import gemmi
    except ImportError as e:
        raise RuntimeError('Install the pinned requirements-structure.txt in the project .venv; see STEP_05B2.md') from e
    I.require(gemmi.__version__ == GEMMI_VERSION, 'Gemmi version differs from tested pin ' + GEMMI_VERSION)
    return gemmi


def file_sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        while chunk := stream.read(1024 * 1024): digest.update(chunk)
    return digest.hexdigest()


def load_input(project):
    path = project / 'config/GSE286387_step05b2_input.json'; b = json.loads(path.read_text())
    I.require(b['status'] == 'STEP05B1_MANIFEST_ACCEPTED_FOR_STEP05B2' and b['protocol_version']=='0.1.0'
              and b['species_taxid']==10090, 'Unaccepted mapping input')
    inputs = [path]
    for name, expected in b['input_sha256'].items():
        f = I.safe_path(project, name)
        I.require(f.is_file() and file_sha(f)==expected, 'Frozen mapping input changed: ' + name); inputs.append(f)
    _, candidates, upstream = D.load_input(project)
    _, sequences, earlier = D.C.load_input(project)
    inputs += upstream + earlier
    ledger = I.read_csv(I.safe_path(project,b['download_manifest_file']))
    archive = {r['path']:r for r in I.read_csv(I.safe_path(project,b['archive_manifest_file']))}
    I.require(len(ledger)==2638 and len({r['file_id'] for r in ledger})==2638 and
              all(r['status']=='downloaded' and r['http_status']=='200' for r in ledger), 'Incomplete acquisition ledger')
    _, expected_links = D.make_plan(candidates)
    links = I.read_csv(I.safe_path(project,b['candidate_links_file']))
    I.require(len(links)==len(candidates)==len(expected_links), 'Candidate link accounting changed')
    for actual, expected in zip(links,expected_links):
        I.require(all(actual[k]==str(v) for k,v in expected.items()),'Candidate links differ from frozen catalogue')
    for r in I.read_csv(I.safe_path(project,b['download_manifest_file']).parent/'input_checksums.csv'):
        f=I.safe_path(project,r['input'])
        I.require(file_sha(f)==r['sha256'],'Acquisition provenance changed: '+r['input']);inputs.append(f)
    return b, sequences, candidates, links, ledger, archive, list(dict.fromkeys(inputs+[Path(__file__).resolve()]))


def verify_raw(project, ledger, archive):
    verified=[]
    for n,r in enumerate(ledger,1):
        for key in ('local_file','metadata_file'):
            name=r[key];f=I.safe_path(project,name)
            I.require(name in archive and f.is_file() and f.stat().st_size==int(archive[name]['bytes']) and
                      file_sha(f)==archive[name]['sha256'], 'Local raw file/sidecar missing or changed: '+name)
        I.require(archive[r['local_file']]['sha256']==r['sha256'] and int(archive[r['local_file']]['bytes'])==int(r['bytes']),
                  'Raw acquisition/archive mismatch')
        meta=json.loads(I.safe_path(project,r['metadata_file']).read_text())
        for key in ('file_id','url','response_url','sha256','retrieved_utc','etag','last_modified'):
            I.require(meta[key]==r[key], 'Raw metadata/acquisition mismatch: '+r['file_id'])
        I.require(meta['http_status']==200 and meta['bytes']==int(r['bytes']), 'Raw metadata byte/status mismatch')
        verified.append(dict(file_id=r['file_id'],kind=r['kind'],local_file=r['local_file'],sha256=r['sha256'],
                             metadata_file=r['metadata_file'],metadata_sha256=archive[r['metadata_file']]['sha256'],status='verified'))
        if n%100==0 or n==len(ledger):print(f'  Raw files verified: {n}/{len(ledger)}',flush=True)
    return verified


def records(block, category):
    table=block.find_mmcif_category(category)
    keys=[str(k).split('.',1)[1] for k in table.tags]
    g=gemmi_module()
    for row in table:
        yield {key:('' if str(value) in ('.','?') else g.cif.as_string(str(value))) for key,value in zip(keys,row)}


def aa(mon):
    letter=gemmi_module().find_tabulated_residue(mon).one_letter_code.upper()
    return letter if letter in 'ACDEFGHIKLMNPQRSTVWYOU' else 'X'


def polymer_index(block):
    entities=defaultdict(dict);bad=set()
    for r in records(block,'_entity_poly_seq.'):
        num=int(r['num']);entity=r['entity_id']
        if num in entities[entity]:bad.add(entity)
        entities[entity][num]=r['mon_id']
    types={r['entity_id']:r['type'] for r in records(block,'_entity_poly.')}
    asym={r['id']:r['entity_id'] for r in records(block,'_struct_asym.')}
    return entities,bad,types,asym


def sequence_map(block,candidate,canonical,label,entity,polymer):
    """Exact sequence or explicit, checked database segments; no heuristic gap alignment."""
    source=candidate['source'];accession=candidate['protein_accession']
    if source=='AlphaFold_DB':
        start,end=int(candidate['sequence_start']),int(candidate['sequence_end'])
        need(sorted(polymer)==list(range(1,end-start+2)), 'PREDICTED_POLYMER_NUMBERING_MISMATCH')
        need(''.join(aa(polymer[n]) for n in sorted(polymer))==canonical[start-1:end], 'PREDICTED_DEPOSITED_SEQUENCE_MISMATCH')
        return {start+n-1:n for n in polymer},'exact_predicted_sequence',[]
    refs={r['id']:r for r in records(block,'_struct_ref.') if r.get('entity_id')==entity and r.get('db_name')=='UNP'}
    exact={k for k,r in refs.items() if r.get('pdbx_db_accession')==accession and not r.get('pdbx_db_isoform')}
    segments=[r for r in records(block,'_struct_ref_seq.') if r.get('ref_id') in exact
              and candidate['chain_id'] in [x.strip() for x in r.get('pdbx_strand_id','').split(',')]]
    mapping={};changes=[]
    diffs=list(records(block,'_struct_ref_seq_dif.'))
    if segments:
        for s in segments:
            need(not any(s.get(k) for k in ('pdbx_seq_align_beg_ins_code','pdbx_seq_align_end_ins_code',
                 'pdbx_db_align_beg_ins_code','pdbx_db_align_end_ins_code')), 'REFERENCE_INSERTION_ALIGNMENT_UNRESOLVED')
            begin,end=int(s['seq_align_beg']),int(s['seq_align_end']);dbbegin,dbend=int(s['db_align_beg']),int(s['db_align_end'])
            need(1<=dbbegin<=dbend<=len(canonical) and begin<=end and end-begin==dbend-dbbegin,
                 'REFERENCE_GAPPED_OR_OUT_OF_RANGE_ALIGNMENT')
            need(all(n in polymer for n in range(begin,end+1)), 'REFERENCE_POLYMER_POSITIONS_MISSING')
            for n in range(begin,end+1):
                c=dbbegin+n-begin;observed=aa(polymer[n]);expected=canonical[c-1]
                if observed!=expected:
                    documented=[r for r in diffs if r.get('align_id')==s['align_id'] and r.get('seq_num')==str(n)
                         and r.get('pdbx_seq_db_seq_num')==str(c) and r.get('pdbx_seq_db_accession_code')==accession
                         and r.get('mon_id') and r.get('db_mon_id') and aa(r['mon_id'])==observed and aa(r['db_mon_id'])==expected]
                    need(len(documented)==1,'REFERENCE_SEQUENCE_DIFFERENCE_UNRESOLVED')
                    changes.append(c)
                need(c not in mapping or mapping[c]==n, 'CONFLICTING_REFERENCE_SEGMENTS')
                mapping[c]=n
        need(len(set(mapping.values()))==len(mapping),'NON_UNIQUE_REFERENCE_MAPPING')
        return mapping,'checked_deposited_UNP_segments',sorted(set(changes))
    need(not refs,'DEPOSITED_UNP_ASSOCIATION_UNRESOLVED')
    need(sorted(polymer)==list(range(1,len(polymer)+1)), 'POLYMER_NUMBERING_NOT_CONTIGUOUS')
    seq=''.join(aa(polymer[n]) for n in sorted(polymer))
    need(bool(seq) and 'X' not in seq,'UNRESOLVED_POLYMER_SEQUENCE')
    first=canonical.find(seq)
    need(first>=0 and canonical.find(seq,first+1)<0,'EXACT_SEQUENCE_PLACEMENT_MISSING_OR_AMBIGUOUS')
    return {first+n:n for n in polymer},'unique_exact_polymer_subsequence',[]


def number(value):
    try:
        n=float(value);return n if math.isfinite(n) else None
    except (ValueError,TypeError):return None


def base_row(candidate,index,pos,file):
    return dict(site_id=f'{candidate["protein_accession"]}:C{pos}',protein_accession=candidate['protein_accession'],
        canonical_cys_position=pos,candidate_row=index,source=candidate['source'],structure_id=candidate['structure_id'],
        candidate_auth_chain=candidate['chain_id'],label_asym_id='',model_id='',label_seq_id='',auth_seq_id='',insertion_code='',
        mapping_method='',mapping_status='held',reason='',deposited_component='',observed_components='',residue_observed=False,
        SG_present=False,SG_alternate_count=0,selected_SG_atom_id='',selected_SG_alt_id='',SG_occupancy='',SG_x='',SG_y='',SG_z='',
        SG_b_iso_or_equiv='',missing_CYS_heavy_atoms='',declared_changed_canonical_positions='',
        structure_file_id=file['file_id'],structure_sha256=file['sha256'],local_file=file['local_file'],
        structural_eligibility='not_assessed',assembly_context='deposited_ASU' if candidate['source']=='PDBe_SIFTS' else 'predicted_model',
        local_quality_status='not_assessed',protocol_version='0.1.0')


def atom_context(block, authors):
    """Index relevant author chains once per file; do not copy unrelated atom payloads."""
    g=gemmi_module();poly=polymer_index(block)
    labels=defaultdict(set);rows=defaultdict(list)
    for r in records(block,'_pdbx_poly_seq_scheme.'):
        for author in [x.strip() for x in r.get('pdb_strand_id','').split(',')]:
            if author in authors:labels[author].add(r['asym_id'])
    table=block.find_mmcif_category('_atom_site.')
    tags={str(t).split('.',1)[1]:i for i,t in enumerate(table.tags)}
    for key in ('label_asym_id','label_seq_id','label_comp_id','label_atom_id','type_symbol','Cartn_x','Cartn_y','Cartn_z'):
        need(key in tags,'ATOM_TABLE_REQUIRED_TAG_MISSING')
    for row in table:
        if str(row[tags['label_seq_id']]) in ('.','?'):continue
        lab=g.cif.as_string(str(row[tags['label_asym_id']]))
        token=str(row[tags['auth_asym_id']]) if 'auth_asym_id' in tags else '.'
        author=lab if token in ('.','?') else g.cif.as_string(token)
        if author in authors:
            labels[author].add(lab);rows[lab].append(row)
    return poly,labels,rows,tags


def map_candidate(block,candidate,index,canonical,file,context=None):
    positions=[n for n,c in enumerate(canonical,1) if c=='C'];out=[]
    auth=candidate['chain_id']
    if context is None:context=atom_context(block,{auth})
    (entities,bad,types,asym),author_labels,atom_rows,tags=context
    labels=author_labels[auth];g=gemmi_module()
    def value(row,key,default=''):
        token=str(row[tags[key]]) if key in tags else default
        return '' if token in ('?','.') else g.cif.as_string(token)
    try:need(bool(labels),'CANDIDATE_AUTHOR_CHAIN_NOT_RESOLVED')
    except MappingHold as e:
        for pos in positions:
            r=base_row(candidate,index,pos,file);r['reason']=str(e);out.append(r)
        return out
    for label in sorted(labels):
        try:
            need(label in asym and asym[label] in entities,'POLYMER_ENTITY_NOT_RESOLVED')
            entity=asym[label];need(entity not in bad,'MICROHETEROGENEOUS_POLYMER_SEQUENCE')
            need(types.get(entity) in ('polypeptide(L)','polypeptide(D)'), 'NON_PROTEIN_OR_UNRESOLVED_POLYMER_TYPE')
            mapping,method,changes=sequence_map(block,candidate,canonical,label,entity,entities[entity])
        except MappingHold as e:
            for pos in positions:
                r=base_row(candidate,index,pos,file);r.update(label_asym_id=label,reason=str(e));out.append(r)
            continue
        target_numbers={mapping[pos] for pos in positions if pos in mapping};atoms=defaultdict(list);models=set()
        for row in atom_rows[label]:
            model=value(row,'pdbx_PDB_model_num') or '1';models.add(model);seq=value(row,'label_seq_id')
            if not seq or int(seq) not in target_numbers:continue
            atom={k:value(row,k) for k in ('id','label_atom_id','label_alt_id','label_comp_id','auth_seq_id',
                  'pdbx_PDB_ins_code','type_symbol','occupancy','Cartn_x','Cartn_y','Cartn_z','B_iso_or_equiv')}
            atoms[(model,int(seq))].append(atom)
        for model in sorted(models or {''}):
            for pos in positions:
                r=base_row(candidate,index,pos,file);r.update(label_asym_id=label,model_id=model,mapping_method=method,
                    declared_changed_canonical_positions=';'.join(map(str,changes)))
                if pos not in mapping:
                    r.update(mapping_status='unassessable',reason='CANONICAL_SITE_OUTSIDE_CHECKED_SEGMENTS');out.append(r);continue
                seq=mapping[pos];found=atoms[(model,seq)];r.update(label_seq_id=seq,deposited_component=entities[entity][seq])
                if not found:
                    r.update(mapping_status='unassessable',reason='RESIDUE_UNOBSERVED');out.append(r);continue
                r['residue_observed']=True;r['observed_components']=';'.join(sorted({a['label_comp_id'] for a in found}))
                identifiers={(a['auth_seq_id'],a['pdbx_PDB_ins_code']) for a in found}
                if len(identifiers)!=1:r['reason']='AUTHOR_RESIDUE_IDENTIFIER_CONFLICT';out.append(r);continue
                r['auth_seq_id'],r['insertion_code']=next(iter(identifiers))
                names={a['label_atom_id'] for a in found};r['missing_CYS_heavy_atoms']=';'.join(sorted({'N','CA','C','O','CB','SG'}-names))
                sg=[a for a in found if a['label_atom_id']=='SG' and a['type_symbol'].upper()=='S']
                r['SG_present']=bool(sg);r['SG_alternate_count']=len(sg)
                if entities[entity][seq]!='CYS' or any(a['label_comp_id']!='CYS' for a in found):
                    r['reason']='TARGET_SUBSTITUTED_OR_MODIFIED';out.append(r);continue
                if not sg:r.update(mapping_status='unassessable',reason='CYS_SG_UNOBSERVED');out.append(r);continue
                if any(number(a['occupancy']) is None or not 0<=number(a['occupancy'])<=1 for a in sg):
                    r['reason']='SG_OCCUPANCY_UNRESOLVED';out.append(r);continue
                ordered=sorted(sg,key=lambda a:(-number(a['occupancy']),a['label_alt_id'],a['id']))
                chosen=ordered[0]
                if number(chosen['occupancy'])<=0:r['reason']='SG_ZERO_OCCUPANCY';out.append(r);continue
                if len({a['label_alt_id'] for a in sg})!=len(sg):r['reason']='DUPLICATE_SG_ALTERNATE_IDENTIFIER';out.append(r);continue
                if any(number(chosen[k]) is None for k in ('Cartn_x','Cartn_y','Cartn_z')):
                    r['reason']='SG_COORDINATES_UNRESOLVED';out.append(r);continue
                r.update(mapping_status='mapped',reason='CYS_SG_MAPPED_QUALITY_PENDING',selected_SG_atom_id=chosen['id'],
                    selected_SG_alt_id=chosen['label_alt_id'],SG_occupancy=number(chosen['occupancy']),
                    SG_x=number(chosen['Cartn_x']),SG_y=number(chosen['Cartn_y']),SG_z=number(chosen['Cartn_z']),
                    SG_b_iso_or_equiv=chosen['B_iso_or_equiv']);out.append(r)
    return out


def site_audit(sequences,mapped,unresolved_metadata=None):
    by_site=defaultdict(list)
    for r in mapped:by_site[r['site_id']].append(r)
    out=[]
    for accession,seq in sorted(sequences.items()):
        for pos,aa_ in enumerate(seq,1):
            if aa_!='C':continue
            sid=f'{accession}:C{pos}';items=by_site[sid];good=sum(r['mapping_status']=='mapped' for r in items)
            pending=any(r['mapping_status']=='held' for r in items) or accession in (unresolved_metadata or set())
            status='held' if good or pending else 'unassessable'
            reason='MAPPED_OPTIONS_AWAIT_LOCAL_QUALITY' if good else 'MAPPING_EVIDENCE_UNRESOLVED' if pending else 'NO_OBSERVED_MAPPED_SG_OPTION' if items else 'NO_METADATA_ELIGIBLE_STRUCTURE'
            out.append(dict(step_id='05b2',object_level='site',object_id=sid,protein_accession=accession,
                canonical_cys_position=pos,status=status,reason_code=reason,mapping_records=len(items),mapped_options=good,
                structural_eligibility='not_assessed',protocol_version='0.1.0'))
    return out


def run(project):
    out=project/'results/05b2_cysteine_mapping';out.mkdir(parents=True,exist_ok=True)
    (out/'SUCCESS.txt').unlink(missing_ok=True)
    stage=Path(tempfile.mkdtemp(prefix='mapping_stage_',dir=out))
    try:
        gemmi=gemmi_module();b,seqs,candidates,links,ledger,archive,inputs=load_input(project)
        print('[1/3] Verifying every local raw file and sidecar...',flush=True)
        verified=verify_raw(project,ledger,archive)
        coordinates={r['file_id']:r for r in ledger if r['kind'].endswith('mmcif')}
        groups=defaultdict(list)
        for index,(c,l) in enumerate(zip(candidates,links),1):
            if c['candidate_status']!='candidate':continue
            ids=[fid for fid in l['file_ids'].split(';') if fid in coordinates]
            I.require(len(ids)==1,'Candidate coordinate link is not unique');groups[ids[0]].append((index,c))
        mapped=[];chain_audit=[]
        print('[2/3] Mapping canonical Cys sites; all candidate models/chains retained...',flush=True)
        for n,(fid,items) in enumerate(sorted(groups.items()),1):
            f=coordinates[fid];block=gemmi.cif.read_file(str(I.safe_path(project,f['local_file']))).sole_block()
            I.require(block.name.casefold()==f['structure_id'].casefold(),'Parsed structure identity differs from acquisition')
            try:context=atom_context(block,{c['chain_id'] for _,c in items});context_error=None
            except MappingHold as e:context=None;context_error=e
            for index,c in items:
                try:
                    if context_error:raise context_error
                    rr=map_candidate(block,c,index,seqs[c['protein_accession']],f,context)
                except MappingHold as e:
                    rr=[]
                    for pos,residue in enumerate(seqs[c['protein_accession']],1):
                        if residue=='C':r=base_row(c,index,pos,f);r['reason']=str(e);rr.append(r)
                mapped.extend(rr)
                chain_audit.append(dict(candidate_row=index,protein_accession=c['protein_accession'],source=c['source'],
                    structure_id=c['structure_id'],auth_chain=c['chain_id'],metadata_status=c['candidate_status'],mapping_records=len(rr),
                    mapped_options=sum(r['mapping_status']=='mapped' for r in rr),
                    mapping_status='mapped_options_quality_pending' if any(r['mapping_status']=='mapped' for r in rr) else 'no_mapped_options',
                    structural_eligibility='not_assessed',protocol_version='0.1.0'))
            del context,block
            if n%25==0 or n==len(groups):print(f'  Coordinate files parsed: {n}/{len(groups)}',flush=True)
        mapped.sort(key=lambda r:(r['site_id'],r['candidate_row'],r['label_asym_id'],r['model_id']))
        for index,c in enumerate(candidates,1):
            if c['candidate_status']=='candidate':continue
            chain_audit.append(dict(candidate_row=index,protein_accession=c['protein_accession'],source=c['source'],
                structure_id=c['structure_id'],auth_chain=c['chain_id'],metadata_status=c['candidate_status'],mapping_records=0,
                mapped_options=0,mapping_status='metadata_'+c['candidate_status'],structural_eligibility='not_assessed',protocol_version='0.1.0'))
        chain_audit.sort(key=lambda r:r['candidate_row'])
        unresolved={c['protein_accession'] for c in candidates if c['candidate_status']=='held'}
        audit=site_audit(seqs,mapped,unresolved);I.require(len(audit)==b['input_sites']==10799,'Canonical site audit incomplete')
        proteins=[]
        for accession in sorted(seqs):
            sites=[r for r in audit if r['protein_accession']==accession]
            proteins.append(dict(protein_accession=accession,canonical_cys_sites=len(sites),
                sites_with_mapped_SG_options=sum(r['mapped_options']>0 for r in sites),
                held_sites=sum(r['status']=='held' for r in sites),unassessable_sites=sum(r['status']=='unassessable' for r in sites),
                structural_eligibility='not_assessed',protocol_version='0.1.0'))
        for name,data in [('site_structure_mapping.csv',mapped),('candidate_mapping_audit.csv',chain_audit),
            ('site_step_audit.csv',audit),('protein_mapping_coverage.csv',proteins),('raw_file_checksums.csv',verified)]:
            I.require(bool(data),'Unexpected empty mapping output');I.write_csv(stage/name,data,list(data[0]))
        inputs=list(dict.fromkeys(inputs+[project/'requirements-structure.txt']))
        I.write_csv(stage/'input_checksums.csv',[dict(input=str(f.relative_to(project)),sha256=file_sha(f)) for f in inputs],['input','sha256'])
        summary=dict(status='STEP05B2_GENERATED_REVIEW_PENDING',implementation_version=VERSION,protocol_version='0.1.0',
            input_proteins=len(seqs),input_sites=len(audit),coordinate_files=len(groups),offered_candidate_rows=len(chain_audit),
            metadata_candidate_chains=sum(r['metadata_status']=='candidate' for r in chain_audit),
            verified_raw_files=len(verified),mapping_records=len(mapped),mapping_record_status_counts=dict(Counter(r['mapping_status'] for r in mapped)),
            mapping_reason_counts=dict(Counter(r['reason'] for r in mapped)),site_status_counts=dict(Counter(r['status'] for r in audit)),
            sites_with_mapped_SG_options=sum(r['mapped_options']>0 for r in audit),proteins_with_mapped_SG_options=sum(r['sites_with_mapped_SG_options']>0 for r in proteins),
            python=platform.python_version(),gemmi=gemmi.__version__,generated_utc=datetime.now(timezone.utc).isoformat())
        (stage/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
        report='\n'.join(['# Step 05b2 — canonical cysteine mapping','','Status: '+summary['status'],'',
            f'- Input proteins: {len(seqs)}; canonical Cys sites: {len(audit)}',
            f'- Verified raw files: {len(verified)}; coordinate files parsed: {len(groups)}',
            f'- Sites with mapped SG options: {summary["sites_with_mapped_SG_options"]}',
            f'- Proteins with mapped SG options: {summary["proteins_with_mapped_SG_options"]}',
            '- Mapping record statuses: '+str(summary['mapping_record_status_counts']),'',
            'Mapping is not structural inclusion or local-quality approval. Every original protein/site remains audited.',
            'All offered metadata-eligible chains, label subchains and deposited models are retained; no best conformation is selected.',
            'Database segments require exact reference identity or explicit checked substitutions; ambiguous/gapped associations are held.',
            'Predicted entity sequence is checked against the frozen canonical interval; author numbering never substitutes for label sequence IDs.',
            'SG selection uses maximum known occupancy, then blank/lexical alternate ID; duplicate alternate IDs are held.',
            'Other atoms are not selected/prepared here. Missing atoms, mutations, modified targets and chemical/assembly context remain explicit.',
            'Next: local confidence/completeness, assembly/mutation context and experimental validation; no SASA/pKa/score calculated.'])+'\n'
        (stage/'mapping_report.md').write_text(report)
        for f in stage.iterdir():shutil.copyfile(f,out/f.name)
        (out/'FAILURE.txt').unlink(missing_ok=True);(out/'SUCCESS.txt').write_text(summary['status']+'\n')
        print('[3/3] Outputs: '+str(out),flush=True);print(report);return summary
    except Exception as e:
        (out/'FAILURE.txt').write_text(str(e)+'\n');raise
    finally:shutil.rmtree(stage)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-dir',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args()
    try:run(args.project_dir.resolve())
    except Exception as e:print('STEP05B2 FAILED: '+str(e),file=sys.stderr);return 1
    return 0


if __name__=='__main__':sys.exit(main())
