#!/usr/bin/env python3
"""Resolve held multi-entry genes using explicit gene-centric representative evidence.

The representative must also be an existing exact-gene UniProt candidate, all
candidate accessions must be accounted for by that gene group, and the original
mouse/canonical/sequence checks must pass. No reviewed/longest-entry fallback.
"""
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import importlib.util
import json
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

SPEC = importlib.util.spec_from_file_location('step03', Path(__file__).with_name('03_map_GSE286387.py'))
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)
API = 'https://www.ebi.ac.uk/proteins/api/genecentric'


def base_accession(value):
    return re.sub(r'-\d+$', '', value)


def stable_gene(value):
    return re.sub(r'\.\d+$', '', value)


def group_links(group):
    """Keep only explicit Ensembl identifiers; names/MGI never substitute for IDs."""
    links = defaultdict(set)
    for record in [group['gene']] + group.get('relatedGene', []):
        gene = stable_gene(record.get('geneName', ''))
        if record.get('geneNameType') == 'Ensembl' and M.GENE.fullmatch(gene):
            links[gene].add(base_accession(record['accession']))
    return links


def choose_representative(gene, options, groups):
    identifiers = {e['primaryAccession'] for e in options}
    relevant = defaultdict(set)
    for group in groups:
        members = group_links(group).get(gene, set())
        if members & identifiers:
            relevant[group['gene']['accession']].update(members)
    if not relevant:
        return None, 'NO_EXPLICIT_GENE_CENTRIC_GROUP', relevant
    if len(relevant) != 1:
        return None, 'MULTIPLE_GENE_CENTRIC_REPRESENTATIVES', relevant
    parent, members = next(iter(relevant.items()))
    if not identifiers <= members:
        return None, 'CANDIDATES_NOT_FULLY_ACCOUNTED_FOR', relevant
    selected = next((e for e in options if e['primaryAccession'] == base_accession(parent)), None)
    if selected is None:
        return None, 'REPRESENTATIVE_LACKS_ORIGINAL_EXACT_GENE_LINK', relevant
    if '-' in parent and parent not in M.displayed_ids(selected):
        return None, 'REPRESENTATIVE_IS_NOT_CURRENT_DISPLAYED_ISOFORM', relevant
    status, reason, entry, isoform, evidence = M.classify(gene, [selected])
    if status != 'pass':
        return None, 'REPRESENTATIVE_' + reason, relevant
    return (entry, isoform, evidence), 'RESOLVED_EXPLICIT_GENE_CENTRIC_REPRESENTATIVE', relevant


def get_page(url, cache):
    u = urllib.parse.urlsplit(url)
    M.require(u.scheme == 'https' and u.netloc == 'www.ebi.ac.uk' and
              u.path == '/proteins/api/genecentric', 'Unexpected gene-centric pagination URL')
    key = M.sha(url.encode())
    body, sidecar = cache / (key + '.json'), cache / (key + '.metadata.json')
    if body.exists() or sidecar.exists():
        M.require(body.exists() and sidecar.exists(), 'Incomplete gene-centric cache pair')
        raw, meta = body.read_bytes(), json.loads(sidecar.read_text())
        M.require(meta['url'] == url and meta['sha256'] == M.sha(raw), 'Corrupt gene-centric cache')
        data = json.loads(raw)
        M.require(isinstance(data, list), 'Invalid cached gene-centric array')
        return data, meta
    for attempt in range(4):
        try:
            request = urllib.request.Request(url, headers={'Accept': 'application/json',
                                                          'User-Agent': 'OIPN-mapping-resolution/1.0.0'})
            with urllib.request.urlopen(request, timeout=90) as response:
                raw = response.read()
                match = re.search(r'<([^>]+)>;\s*rel="next"', response.headers.get('Link', ''))
                meta = dict(url=url, sha256=M.sha(raw), retrieved_utc=datetime.now(timezone.utc).isoformat(),
                            total_results=int(response.headers['X-Pagination-TotalRecords']),
                            next_url=match.group(1) if match else '',
                            reference_release=response.headers.get('X-UniProt-Release', 'not_supplied_by_service'))
            data = json.loads(raw)
            M.require(isinstance(data, list) and all('gene' in g for g in data), 'Invalid gene-centric response')
            body.write_bytes(raw); sidecar.write_text(json.dumps(meta, indent=2) + '\n')
            return data, meta
        except urllib.error.HTTPError as error:
            if error.code not in (429, 500, 502, 503, 504) or attempt == 3:
                raise
            retry = error.headers.get('Retry-After', '')
            time.sleep(min(30, int(retry)) if retry.isdigit() else 2 ** attempt)
        except (urllib.error.URLError, TimeoutError):
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    raise RuntimeError('Gene-centric request retries exhausted')


def fetch_batch(accessions, cache, index):
    url = API + '?' + urllib.parse.urlencode({'offset': 0, 'size': 100, 'accession': ','.join(accessions)})
    groups, sources, seen, totals = [], [], set(), set()
    while url:
        M.require(url not in seen, 'Repeated gene-centric pagination URL'); seen.add(url)
        data, meta = get_page(url, cache)
        groups.extend(data); totals.add(meta['total_results'])
        sources.append(dict(meta, batch_id=index, page_records=len(data)))
        url = meta['next_url']
    M.require(len(totals) == 1 and len(groups) == next(iter(totals)), 'Incomplete gene-centric pagination')
    return groups, sources


def load_baseline(project, binding, genes):
    directory = project / 'results/03_protein_mapping'
    M.require((directory / 'SUCCESS.txt').is_file() and not (directory / 'FAILURE.txt').exists(),
              'Run Step 03 successfully before Step 03b')
    M.require((directory / 'SUCCESS.txt').read_text().strip() == 'STEP03_GENERATED_REVIEW_PENDING',
              'Unexpected Step 03 success marker')
    for record in M.read_csv(directory / 'input_checksums.csv'):
        path = project / record['input']
        M.require(path.is_file() and M.sha(path.read_bytes()) == record['sha256'],
                  'Step 03 input or helper-script checksum changed: ' + record['input'])
    rows = M.read_csv(directory / 'gene_protein_mapping.csv')
    M.require(len(rows) == len(genes) and {r['stable_gene_id'] for r in rows} ==
              {r['stable_gene_id'] for r in genes}, 'Step 03 gene universe mismatch')
    cache = project / 'data/raw/UniProt/step03' / binding['discovery_sha256']
    sources = M.read_csv(directory / 'source_manifest.csv')
    entries, releases = {}, set()
    for source in sources:
        key = M.sha(source['url'].encode())
        raw = (cache / (key + '.json')).read_bytes()
        meta = json.loads((cache / (key + '.metadata.json')).read_text())
        page, meta = M.validated_page(raw, meta, source['url'])
        M.require(M.sha(raw) == source['sha256'] and meta['release'] == source['release'], 'Step 03 source provenance differs')
        releases.add(meta['release'])
        for entry in page['results']:
            accession = entry['primaryAccession']
            M.require(accession not in entries or entry == entries[accession], 'Baseline entry disagreement')
            entries[accession] = entry
    M.require(len(releases) == 1, 'Baseline mixed releases')
    by_gene = defaultdict(list)
    for entry in entries.values():
        for gene in M.gene_refs(entry):
            by_gene[gene].append(entry)
    for row in rows:
        options = sorted(by_gene[row['stable_gene_id']], key=lambda e: e['primaryAccession'])
        status, reason, entry, isoform, evidence = M.classify(row['stable_gene_id'], options)
        M.require((row['mapping_status'], row['reason']) == (status, reason) and
                  row['candidate_accessions'] == ';'.join(e['primaryAccession'] for e in options),
                  'Baseline mapping status/candidates changed')
        if entry:
            M.require(row['protein_accession'] == entry['primaryAccession'] and
                      row['sequence_checksum'] == M.sha(entry['sequence']['value'].encode()), 'Baseline sequence mismatch')
        else:
            M.require(not row['protein_accession'] and not row['sequence_checksum'], 'Unaccepted baseline has a claimed sequence')
    return rows, by_gene, next(iter(releases)), directory


def update_row(row, selection, reason):
    entry, isoform, evidence = selection
    gene = row['stable_gene_id']; refs = M.gene_refs(entry)[gene]; seq = entry['sequence']['value']
    row.update(mapping_status='pass', reason=reason, protein_accession=entry['primaryAccession'],
               isoform_id=isoform, canonical_status='UniProtKB_displayed_sequence', species_taxid=10090,
               sequence_length=len(seq), sequence_checksum=M.sha(seq.encode()),
               canonical_gene_association=evidence, entry_type=entry.get('entryType', ''),
               entry_version=entry.get('entryAudit', {}).get('entryVersion', ''),
               sequence_version=entry.get('entryAudit', {}).get('sequenceVersion', ''),
               ensembl_gene_ids_versioned=';'.join(sorted({r['gene_id_versioned'] for r in refs})),
               ensembl_transcript_ids=';'.join(sorted({r['transcript_id'] for r in refs})),
               ensembl_protein_ids=';'.join(sorted({r['protein_id'] for r in refs})))


def run(project):
    out = project / 'results/03b_mapping_resolution'; out.mkdir(parents=True, exist_ok=True)
    (out / 'SUCCESS.txt').unlink(missing_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='resolution_stage_', dir=out))
    try:
        binding, genes, inputs = M.load_input(project)
        rows, by_gene, release, baseline = load_baseline(project, binding, genes)
        held = [r for r in rows if r['reason'] == 'MULTIPLE_MOUSE_UNIPROT_ENTRIES']
        original_pass = sum(r['mapping_status'] == 'pass' for r in rows)
        accessions = sorted({e['primaryAccession'] for r in held for e in by_gene[r['stable_gene_id']]})
        cache = project / 'data/raw/gene_centric/step03b' / binding['discovery_sha256']
        cache.mkdir(parents=True, exist_ok=True)
        batches = [accessions[i:i + 35] for i in range(0, len(accessions), 35)]
        groups, sources = [], []
        print(f'[1/4] Checking {len(held)} held genes using {len(batches)} gene-centric batches...', flush=True)
        with ThreadPoolExecutor(max_workers=3) as pool:
            jobs = [pool.submit(fetch_batch, batch, cache, i + 1) for i, batch in enumerate(batches)]
            for done, job in enumerate(as_completed(jobs), 1):
                data, provenance = job.result(); groups.extend(data); sources.extend(provenance)
                print(f'  Gene-centric batches complete: {done}/{len(jobs)}', flush=True)
        known_releases = {s['reference_release'] for s in sources if s['reference_release'] != 'not_supplied_by_service'}
        M.require(not known_releases or known_releases == {release},
                  'Gene-centric declared release differs from the frozen UniProt sequence release')
        # Repeated groups across accession queries are normal; deduplicate identical evidence.
        groups = list({json.dumps(g, sort_keys=True): g for g in groups}.values())
        evidence_rows, audit, sequences = [], [], {}
        print('[2/4] Requiring one explicit representative and complete candidate accounting...', flush=True)
        for row in rows:
            gene = row['stable_gene_id']; row['step03_original_status'] = row['mapping_status']
            row['resolution_method'] = 'baseline_unchanged'; row['resolution_version'] = '1.0.0'
            if row['reason'] == 'MULTIPLE_MOUSE_UNIPROT_ENTRIES':
                options = by_gene[gene]
                selection, reason, relevant = choose_representative(gene, options, groups)
                for parent, members in sorted(relevant.items()):
                    evidence_rows.append(dict(stable_gene_id=gene, representative_identifier=parent,
                                              member_accessions=';'.join(sorted(members)),
                                              covers_all_original_candidates={e['primaryAccession'] for e in options} <= members,
                                              resolution_reason=reason))
                row['resolution_method'] = 'explicit_gene_centric_evidence'
                row['reason'] = reason
                if selection:
                    update_row(row, selection, reason)
            if row['mapping_status'] == 'pass':
                entry = next(e for e in by_gene[gene] if e['primaryAccession'] == row['protein_accession'])
                seq = entry['sequence']['value']
                M.require(M.sha(seq.encode()) == row['sequence_checksum'], 'Final sequence hash mismatch')
                sequences[entry['primaryAccession']] = seq
            audit.append(dict(step_id='03', object_level='gene', object_id=gene, status=row['mapping_status'],
                              reason_code=row['reason'], reason_text=row['reason'], input_source=binding['discovery_file'],
                              input_checksum=binding['discovery_sha256'], protocol_version=binding['protocol_version']))
        print('[3/4] Saving combined mapping, FASTA and evidence; baseline preserved...', flush=True)
        M.write_csv(stage / 'gene_protein_mapping.csv', rows, list(rows[0]))
        M.write_csv(stage / 'step_audit.csv', audit, list(audit[0]))
        M.write_csv(stage / 'gene_centric_evidence.csv', evidence_rows,
                    ['stable_gene_id', 'representative_identifier', 'member_accessions', 'covers_all_original_candidates', 'resolution_reason'])
        with (stage / 'canonical_sequences.fasta').open('w') as handle:
            for accession, seq in sorted(sequences.items()):
                associated = sorted(r['stable_gene_id'] for r in rows if r['protein_accession'] == accession)
                handle.write('>' + accession + ' taxid=10090 genes=' + ','.join(associated) + '\n')
                handle.write('\n'.join(seq[i:i + 60] for i in range(0, len(seq), 60)) + '\n')
        sources.sort(key=lambda r: (r['batch_id'], r['url']))
        M.write_csv(stage / 'source_manifest.csv', sources,
                    ['url', 'sha256', 'retrieved_utc', 'total_results', 'next_url', 'reference_release', 'batch_id', 'page_records'])
        inputs += [baseline / 'gene_protein_mapping.csv', baseline / 'source_manifest.csv', Path(__file__).resolve()]
        M.write_csv(stage / 'input_checksums.csv', [{'input': str(p.relative_to(project) if project in p.parents else p),
                                                   'sha256': M.sha(p.read_bytes())} for p in inputs], ['input', 'sha256'])
        passed = sum(r['mapping_status'] == 'pass' for r in rows)
        summary = dict(status='STEP03B_GENERATED_REVIEW_PENDING', input_genes=len(rows), baseline_passed_genes=original_pass,
                       resolved_additional_genes=passed - original_pass, passed_genes=passed,
                       unique_accepted_proteins=len(sequences), mapped_gene_fraction=passed / len(rows),
                       status_counts=dict(Counter(r['mapping_status'] for r in rows)),
                       reason_counts=dict(Counter(r['reason'] for r in rows)), uniprot_sequence_release=release,
                       gene_centric_release_headers=sorted({r['reference_release'] for r in sources}),
                       protocol_version=binding['protocol_version'], resolution_version='1.0.0',
                       discovery_sha256=binding['discovery_sha256'], python=platform.python_version())
        (stage / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        report = ['# Step 03b — evidence-based mapping resolution', '', 'Status: STEP03B_GENERATED_REVIEW_PENDING', '',
                  f'- Input genes: {len(rows)}', f'- Baseline passed: {original_pass}',
                  f'- Additional resolved: {passed - original_pass}', f'- Total passed genes: {passed} ({passed/len(rows):.1%})',
                  f'- Unique accepted proteins: {len(sequences)}', f'- Status counts: {summary["status_counts"]}',
                  f'- UniProt sequence release: {release}', '',
                  'Gene-centric representatives require explicit Ensembl gene links, complete candidate accounting, and existing exact-gene mouse UniProt sequence evidence.',
                  'No entry is chosen merely because it is reviewed, longer or biologically interesting. Step 03 baseline is preserved.',
                  'Gene-centric service may not provide a release header; its raw response SHA256 and retrieval time are recorded separately.',
                  'The six missing-current-link genes and any fragment remain audited; this substep does not infer historical-ID replacements or repair incomplete sequences.',
                  'Coverage is descriptive, not a threshold to force mapping. Workstation reproduction, evidence and unresolved cases require review before Step 04.']
        (stage / 'resolution_report.md').write_text('\n'.join(report) + '\n')
        for path in stage.iterdir():
            shutil.copyfile(path, out / path.name)
        (out / 'FAILURE.txt').unlink(missing_ok=True)
        (out / 'SUCCESS.txt').write_text('STEP03B_GENERATED_REVIEW_PENDING\n')
        print('[4/4] Outputs: ' + str(out), flush=True); print('\n'.join(report[:11]))
        return summary
    except Exception as error:
        (out / 'FAILURE.txt').write_text(str(error) + '\n'); raise
    finally:
        shutil.rmtree(stage)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-dir', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        run(args.project_dir.resolve())
    except Exception as error:
        print('STEP03B FAILED: ' + str(error), file=sys.stderr); return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
