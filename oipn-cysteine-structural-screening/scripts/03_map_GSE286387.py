#!/usr/bin/env python3
"""Checksum-bound mouse gene -> UniProtKB displayed/canonical sequence audit.

Python standard library only. Current exact Ensembl GeneId cross-references are
used; historical-ID rescue, symbols and arbitrary multi-entry selection are not.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import platform
import re
import shutil
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

API = 'https://rest.uniprot.org/uniprotkb/search'
AA = set('ACDEFGHIKLMNPQRSTVWYOU')
GENE = re.compile(r'^ENSMUSG\d{11}$')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read_csv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def gene_refs(entry):
    result = defaultdict(list)
    for ref in entry.get('uniProtKBCrossReferences', []):
        if ref.get('database') != 'Ensembl':
            continue
        props = {p['key']: p['value'] for p in ref.get('properties', [])}
        raw = props.get('GeneId', '')
        stable = re.sub(r'\.\d+$', '', raw)
        if GENE.fullmatch(stable):
            result[stable].append({'gene_id_versioned': raw, 'transcript_id': ref['id'],
                                   'protein_id': props.get('ProteinId', ''),
                                   'isoform_id': ref.get('isoformId', '')})
    return result


def displayed_ids(entry):
    return sorted({ident for comment in entry.get('comments', [])
                   if comment.get('commentType') == 'ALTERNATIVE PRODUCTS'
                   for iso in comment.get('isoforms', [])
                   if iso.get('isoformSequenceStatus') == 'Displayed'
                   for ident in iso.get('isoformIds', [])})


def sequence_quality(entry):
    seq = entry.get('sequence', {}).get('value', '')
    if not seq:
        return 'MISSING_SEQUENCE'
    if set(seq) - AA:
        return 'AMBIGUOUS_OR_INVALID_SEQUENCE'
    if entry.get('sequence', {}).get('length') != len(seq):
        return 'SEQUENCE_LENGTH_MISMATCH'
    md5 = entry.get('sequence', {}).get('md5')
    if md5 and hashlib.md5(seq.encode()).hexdigest().upper() != md5.upper():
        return 'SEQUENCE_MD5_MISMATCH'
    if 'fragment' in entry.get('proteinDescription', {}).get('flag', '').lower():
        return 'FRAGMENT_SEQUENCE'
    return ''


def canonical_evidence(entry, refs):
    ids = displayed_ids(entry)
    if len(ids) > 1:
        return '', '', 'MULTIPLE_DISPLAYED_ISOFORMS'
    if any(r['isoform_id'] and r['isoform_id'] in ids for r in refs):
        return ids[0], 'explicit_displayed_isoform_cross_reference', ''
    if any(not r['isoform_id'] for r in refs):
        return ids[0] if ids else '', 'entry_level_cross_reference_isoform_unspecified', ''
    return '', '', 'ONLY_NONCANONICAL_OR_UNRESOLVED_ISOFORM_XREF'


def classify(gene, candidates):
    """Return one terminal gene status; do not resolve multiple accessions by quality."""
    if not candidates:
        return 'held', 'NO_CURRENT_EXACT_ENSEMBL_XREF', None, '', ''
    mouse = [e for e in candidates if e.get('organism', {}).get('taxonId') == 10090]
    if not mouse:
        return 'held', 'WRONG_SPECIES', None, '', ''
    if len(mouse) != 1:
        return 'held', 'MULTIPLE_MOUSE_UNIPROT_ENTRIES', None, '', ''
    entry = mouse[0]
    require(gene in gene_refs(entry), 'Candidate lacks exact source gene cross-reference')
    accession = entry.get('primaryAccession', '')
    if not re.fullmatch(r'[A-Z0-9]{6,10}', accession):
        return 'held', 'UNRESOLVED_PRIMARY_ACCESSION', None, '', ''
    isoform, evidence, reason = canonical_evidence(entry, gene_refs(entry)[gene])
    if reason:
        return 'held', reason, None, '', ''
    reason = sequence_quality(entry)
    if reason:
        return 'unassessable', reason, None, '', ''
    return 'pass', 'PASS_CURRENT_EXACT_XREF_DISPLAYED_SEQUENCE', entry, isoform, evidence


def validated_page(raw, meta, url):
    require(meta.get('url') == url and meta.get('sha256') == sha(raw), 'Corrupt/mismatched UniProt cache')
    require(meta.get('release') and meta.get('release_date'), 'Missing UniProt release headers')
    page = json.loads(raw)
    require(isinstance(page.get('results'), list), 'UniProt response lacks results array')
    return page, meta


def get_page(url, cache):
    parsed = urllib.parse.urlsplit(url)
    require(parsed.scheme == 'https' and parsed.netloc == 'rest.uniprot.org' and
            parsed.path == '/uniprotkb/search', 'Unexpected pagination URL')
    key = sha(url.encode())
    data_path, meta_path = cache / (key + '.json'), cache / (key + '.metadata.json')
    if data_path.exists() or meta_path.exists():
        require(data_path.exists() and meta_path.exists(), 'Incomplete cache pair; remove this pair before retrying')
        return validated_page(data_path.read_bytes(), json.loads(meta_path.read_text()), url)
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'OIPN-structural-screening/0.1.0',
                                                          'Accept': 'application/json'})
            with urllib.request.urlopen(request, timeout=90) as response:
                raw = response.read()
                link = response.headers.get('Link', '')
                match = re.search(r'<([^>]+)>;\s*rel="next"', link)
                meta = {'url': url, 'sha256': sha(raw), 'retrieved_utc': datetime.now(timezone.utc).isoformat(),
                        'release': response.headers.get('X-UniProt-Release'),
                        'release_date': response.headers.get('X-UniProt-Release-Date'),
                        'total_results': int(response.headers['X-Total-Results']),
                        'next_url': match.group(1) if match else ''}
            page, meta = validated_page(raw, meta, url)
            # Completed HTTP/schema checks precede publication of a cached response.
            data_path.write_bytes(raw)
            meta_path.write_text(json.dumps(meta, indent=2) + '\n')
            return page, meta
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504) or attempt == 3:
                raise
            retry = error.headers.get('Retry-After', '')
            time.sleep(min(30, int(retry)) if retry.isdigit() else 2 ** attempt)
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError('UniProt request retries exhausted')


def fetch_batch(genes, cache, batch_id):
    query = 'organism_id:10090 AND (' + ' OR '.join('xref:ensembl-' + g for g in genes) + ')'
    url = API + '?' + urllib.parse.urlencode({'query': query, 'format': 'json', 'size': 500})
    entries, provenance, seen, totals = [], [], set(), set()
    while url:
        require(url not in seen, 'Repeated UniProt pagination URL')
        seen.add(url)
        page, meta = get_page(url, cache)
        totals.add(meta['total_results'])
        entries.extend(page['results'])
        provenance.append(dict(meta, batch_id=batch_id, page_entries=len(page['results'])))
        url = meta['next_url']
    require(len(totals) == 1 and len(entries) == next(iter(totals)), 'Incomplete UniProt pagination')
    accessions = [e.get('primaryAccession') for e in entries]
    require(None not in accessions and len(set(accessions)) == len(accessions), 'Duplicate/missing page accessions')
    return entries, provenance


def load_input(project):
    binding_path = project / 'config/GSE286387_step03_input.json'
    binding = json.loads(binding_path.read_text())
    cfg_path = project / 'config/analysis_parameters.json'
    cfg = json.loads(cfg_path.read_text())
    require(binding['status'] == 'STEP02_SNAPSHOT_ACCEPTED_FOR_STEP03' and
            binding['species_taxid'] == cfg['scope']['taxonomy_id'] == 10090 and
            binding['protocol_version'] == cfg['protocol_version'] == '0.1.0', 'Unaccepted or mismatched Step 03 input')
    require(cfg['sequence']['species_must_match'] is True and
            cfg['sequence']['primary_isoform'] == 'UniProt_canonical' and
            cfg['sequence']['ambiguous_mapping'] == 'hold_for_review', 'Unsupported sequence mapping policy')
    input_path = (project / binding['discovery_file']).resolve()
    require(project in input_path.parents, 'Discovery path escapes project')
    require(sha(input_path.read_bytes()) == binding['discovery_sha256'], 'Accepted discovery CSV changed')
    genes = read_csv(input_path)
    ids = [r['stable_gene_id'] for r in genes]
    require(len(ids) == binding['n_unique_genes'] and len(set(ids)) == len(ids) and
            all(GENE.fullmatch(x) for x in ids), 'Invalid/duplicate discovery gene IDs')
    require(all(r['biotype'] == 'protein_coding' and r['passes_to_step03'] == 'TRUE' and
                math.isfinite(float(r['padj'])) and 0 <= float(r['padj']) < .05 for r in genes),
            'Discovery violates registered coding/FDR selection')
    require(Counter(r['direction'] for r in genes) == {'up': binding['n_up'], 'down': binding['n_down']},
            'Discovery direction counts differ')
    return binding, genes, [binding_path, cfg_path, input_path, Path(__file__).resolve()]


def run(project):
    out = project / 'results/03_protein_mapping'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'SUCCESS.txt').unlink(missing_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='mapping_stage_', dir=out))
    try:
        binding, genes, inputs = load_input(project)
        cache = project / 'data/raw/UniProt/step03' / binding['discovery_sha256']
        cache.mkdir(parents=True, exist_ok=True)
        ids = sorted(r['stable_gene_id'] for r in genes)
        batches = [ids[i:i + 25] for i in range(0, len(ids), 25)]
        entries, provenance = {}, []
        print(f'[1/4] Accepted {len(ids)} genes; retrieving {len(batches)} cached/paginated UniProt batches...', flush=True)
        with ThreadPoolExecutor(max_workers=3) as pool:
            jobs = [pool.submit(fetch_batch, group, cache, i + 1) for i, group in enumerate(batches)]
            for done, job in enumerate(as_completed(jobs), 1):
                batch_entries, sources = job.result()
                provenance.extend(sources)
                for entry in batch_entries:
                    accession = entry['primaryAccession']
                    require(accession not in entries or entries[accession] == entry,
                            'Same accession changed between queries; mixed reference snapshot')
                    entries[accession] = entry
                print(f'  UniProt batches complete: {done}/{len(jobs)}', flush=True)
        releases = {(p['release'], p['release_date']) for p in provenance}
        require(len(releases) == 1, 'Mixed UniProt releases; archive cache and rerun in one release')
        release, release_date = next(iter(releases))
        by_gene = defaultdict(list)
        for entry in entries.values():
            for gene in gene_refs(entry):
                if gene in ids:
                    by_gene[gene].append(entry)
        print('[2/4] Auditing exact gene associations and canonical sequence eligibility...', flush=True)
        mapping, candidates, audit, accepted = [], [], [], {}
        for gene_row in genes:
            gene = gene_row['stable_gene_id']
            options = sorted(by_gene[gene], key=lambda e: e['primaryAccession'])
            status, reason, entry, isoform, evidence = classify(gene, options)
            row = dict(stable_gene_id=gene, symbol=gene_row['symbol'], protein_accession='', isoform_id='',
                       canonical_status='', species_taxid='', sequence_length='', sequence_checksum='',
                       reference_release=release, reference_release_date=release_date,
                       input_annotation_reference=binding['annotation_reference'], mapping_status=status,
                       reason=reason, n_candidate_accessions=len(options), candidate_accessions=';'.join(e['primaryAccession'] for e in options),
                       canonical_gene_association='', entry_type='', entry_version='', sequence_version='',
                       ensembl_gene_ids_versioned='', ensembl_transcript_ids='', ensembl_protein_ids='',
                       protocol_version=binding['protocol_version'])
            for option in options:
                refs = gene_refs(option)[gene]
                candidate_iso, candidate_evidence, candidate_reason = canonical_evidence(option, refs)
                candidates.append(dict(stable_gene_id=gene, protein_accession=option['primaryAccession'],
                                       species_taxid=option.get('organism', {}).get('taxonId', ''),
                                       entry_type=option.get('entryType', ''), displayed_isoform_id=candidate_iso,
                                       canonical_gene_association=candidate_evidence, canonical_xref_issue=candidate_reason,
                                       sequence_issue=sequence_quality(option),
                                       ensembl_cross_references_json=json.dumps(refs, sort_keys=True),
                                       reference_release=release, protocol_version=binding['protocol_version']))
            if entry:
                seq = entry['sequence']['value']
                refs = gene_refs(entry)[gene]
                row.update(protein_accession=entry['primaryAccession'], isoform_id=isoform,
                           canonical_status='UniProtKB_displayed_sequence', species_taxid=10090,
                           sequence_length=len(seq), sequence_checksum=sha(seq.encode()),
                           canonical_gene_association=evidence, entry_type=entry.get('entryType', ''),
                           entry_version=entry.get('entryAudit', {}).get('entryVersion', ''),
                           sequence_version=entry.get('entryAudit', {}).get('sequenceVersion', ''),
                           ensembl_gene_ids_versioned=';'.join(sorted({r['gene_id_versioned'] for r in refs})),
                           ensembl_transcript_ids=';'.join(sorted({r['transcript_id'] for r in refs})),
                           ensembl_protein_ids=';'.join(sorted({r['protein_id'] for r in refs})))
                require(entry['primaryAccession'] not in accepted or accepted[entry['primaryAccession']] == seq,
                        'Same accepted accession has different sequences')
                accepted[entry['primaryAccession']] = seq
            mapping.append(row)
            audit.append(dict(step_id='03', object_level='gene', object_id=gene, status=status,
                              reason_code=reason, reason_text=reason, input_source=binding['discovery_file'],
                              input_checksum=binding['discovery_sha256'], protocol_version=binding['protocol_version']))
        print('[3/4] Saving complete gene/candidate audits and unique-protein FASTA...', flush=True)
        write_csv(stage / 'gene_protein_mapping.csv', mapping, list(mapping[0]))
        candidate_fields = ['stable_gene_id', 'protein_accession', 'species_taxid', 'entry_type', 'displayed_isoform_id',
                            'canonical_gene_association', 'canonical_xref_issue', 'sequence_issue',
                            'ensembl_cross_references_json', 'reference_release', 'protocol_version']
        write_csv(stage / 'mapping_candidates.csv', candidates, candidate_fields)
        write_csv(stage / 'step_audit.csv', audit, list(audit[0]))
        with (stage / 'canonical_sequences.fasta').open('w') as handle:
            for accession, seq in sorted(accepted.items()):
                associated = sorted(r['stable_gene_id'] for r in mapping if r['protein_accession'] == accession)
                handle.write('>' + accession + ' taxid=10090 genes=' + ','.join(associated) + '\n')
                handle.write('\n'.join(seq[i:i + 60] for i in range(0, len(seq), 60)) + '\n')
        provenance.sort(key=lambda r: (r['batch_id'], r['url']))
        write_csv(stage / 'source_manifest.csv', provenance, list(provenance[0]))
        write_csv(stage / 'input_checksums.csv', [{'input': str(p.relative_to(project) if project in p.parents else p),
                                                'sha256': sha(p.read_bytes())} for p in inputs], ['input', 'sha256'])
        summary = dict(status='STEP03_GENERATED_REVIEW_PENDING', input_genes=len(genes),
                       mapping_status_counts=dict(Counter(r['mapping_status'] for r in mapping)),
                       reason_counts=dict(Counter(r['reason'] for r in mapping)),
                       passed_genes=sum(r['mapping_status'] == 'pass' for r in mapping),
                       unique_accepted_proteins=len(accepted), candidate_accessions=len(entries),
                       reference_release=release, reference_release_date=release_date,
                       discovery_sha256=binding['discovery_sha256'], protocol_version=binding['protocol_version'],
                       python=platform.python_version(), generated_utc=datetime.now(timezone.utc).isoformat())
        (stage / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        report = ['# Step 03 — mouse canonical protein mapping', '', 'Status: STEP03_GENERATED_REVIEW_PENDING', '',
                  f'- Input genes: {len(genes)}', f'- Passed genes: {summary["passed_genes"]}',
                  f'- Unique accepted proteins: {len(accepted)}', f'- Status counts: {summary["mapping_status_counts"]}',
                  f'- UniProt release: {release} ({release_date})', '',
                  'Exact current Ensembl GeneId cross-references only; historical release 76 IDs are not silently rescued by symbols.',
                  'Multiple mouse accessions are held, including reviewed/unreviewed alternatives; none is chosen by length, score or review status.',
                  'Canonical means UniProtKB displayed sequence; entry-level unspecified isoform association is explicitly distinguished from an explicit displayed-isoform cross-reference.',
                  'RNA does not establish the expressed isoform. Canonical sequences may be precursors; no processing or Cys filtering occurs here.',
                  'Missing/multiple mappings are held, not evidence of absent protein or oxidative resistance.',
                  'Review held/unassessable reasons and full gene/candidate audit before Step 04.']
        (stage / 'mapping_report.md').write_text('\n'.join(report) + '\n')
        for path in stage.iterdir():
            shutil.copyfile(path, out / path.name)
        (out / 'FAILURE.txt').unlink(missing_ok=True)
        (out / 'SUCCESS.txt').write_text('STEP03_GENERATED_REVIEW_PENDING\n')
        print('[4/4] Outputs: ' + str(out), flush=True)
        print('\n'.join(report[:10]))
        return summary
    except Exception as error:
        (out / 'FAILURE.txt').write_text(str(error) + '\n')
        raise
    finally:
        shutil.rmtree(stage)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-dir', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        run(args.project_dir.resolve())
    except Exception as error:
        print('STEP03 FAILED: ' + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
