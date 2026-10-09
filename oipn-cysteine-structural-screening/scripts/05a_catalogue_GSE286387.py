#!/usr/bin/env python3
"""Catalogue structural candidates; coordinates/local quality are Step 05b work."""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import importlib.util
import io
import json
import math
from pathlib import Path
import platform
import re
import shutil
import sys
import tempfile
import tarfile
from datetime import datetime, timezone
import time
import urllib.error
import urllib.parse
import urllib.request

SPEC = importlib.util.spec_from_file_location('inventory', Path(__file__).with_name('04_inventory_GSE286387.py'))
I = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(I)
VERSION = '1.0.2'
PDBE = 'https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/'
AFDB = 'https://alphafold.ebi.ac.uk/api/prediction/'


class CatalogueHTTPError(urllib.error.HTTPError):
    """Retain the requested endpoint and a bounded access-error diagnostic."""

    def __init__(self, requested_url, original):
        raw = original.read(1025)
        excerpt = raw[:1024]
        headers = original.headers or {}
        super().__init__(original.geturl(), original.code, original.reason,
                         original.headers, io.BytesIO(excerpt))
        self.details = dict(
            status='HTTP_REQUEST_FAILED', requested_url=requested_url,
            response_url=original.geturl(), http_status=original.code,
            reason=str(original.reason)[:200],
            response_headers={key: str(headers.get(key, ''))[:256] for key in
                              ('Server', 'Content-Type', 'Retry-After', 'Via', 'X-Mitmproxy-Blocked-Reason')},
            response_excerpt=excerpt.decode('utf-8', errors='replace'),
            excerpt_truncated=len(raw) > 1024,
            retrieved_utc=datetime.now(timezone.utc).isoformat())

    def __str__(self):
        return f'HTTP {self.code}: {self.reason}; requested_url={self.details["requested_url"]}'


def load_input(project):
    path = project / 'config/GSE286387_step05_input.json'
    binding = json.loads(path.read_text())
    I.require(binding['status'] == 'STEP04_SNAPSHOT_ACCEPTED_FOR_STEP05' and
              binding['protocol_version'] == '0.1.0' and binding['species_taxid'] == 10090,
              'Step 04 review/species not accepted')
    inputs = [path]
    for relative, expected in binding['input_sha256'].items():
        f = I.safe_path(project, relative)
        I.require(f.is_file() and I.sha(f.read_bytes()) == expected, 'Frozen Step 05 input changed: ' + relative)
        inputs.append(f)
    for relative in (binding['fasta_file'], binding['site_file'], binding['protein_file']):
        I.require(relative in binding['input_sha256'], 'Primary input is not checksum-bound')
    old_binding, rows, originals, associations, upstream = I.load_input(project)
    sequences, metadata = I.read_fasta(I.safe_path(project, binding['fasta_file']))
    proteins = I.read_csv(I.safe_path(project, binding['protein_file']))
    positive = {r['protein_accession'] for r in proteins if r['inventory_status'] == 'pass'}
    I.require(len(sequences) == binding['input_proteins'] and set(sequences) == positive,
              'Cys-positive protein universe mismatch')
    for accession, seq in sequences.items():
        I.require(seq == originals[accession] and metadata[accession] == associations[accession],
                  'Canonical sequence/gene association changed')
    sites = I.read_csv(I.safe_path(project, binding['site_file']))
    expected = {(a, n) for a, seq in sequences.items() for n, residue in enumerate(seq, 1) if residue == 'C'}
    actual = {(r['protein_accession'], int(r['canonical_cys_position'])) for r in sites}
    I.require(len(sites) == len(actual) == binding['input_sites'] and actual == expected,
              'Canonical Cys site universe mismatch')
    I.require(all(r['site_id'] == f'{r["protein_accession"]}:C{r["canonical_cys_position"]}' for r in sites),
              'Cys site identifier mismatch')
    return binding, sequences, inputs + upstream


def fetch(url, cache, offline=False):
    """Preserve 200/404 response bytes. Network/server errors are never absence."""
    I.require(url.startswith(PDBE) or url.startswith(AFDB), 'Unexpected catalogue API URL')
    key = I.sha(url.encode())
    body, sidecar = cache / (key + '.response'), cache / (key + '.metadata.json')
    if body.exists() or sidecar.exists():
        I.require(body.exists() and sidecar.exists(), 'Incomplete structural metadata cache pair')
        raw, meta = body.read_bytes(), json.loads(sidecar.read_text())
        I.require(meta['url'] == url and meta['sha256'] == I.sha(raw) and meta['http_status'] in (200, 404),
                  'Corrupt structural metadata cache')
    else:
        I.require(not offline, 'Offline cache missing: ' + url)
        for attempt in range(4):
            try:
                request = urllib.request.Request(url, headers={'Accept': 'application/json',
                                                              'User-Agent': 'OIPN-structure-catalogue/1.0.0'})
                with urllib.request.urlopen(request, timeout=60) as response:
                    raw, code, headers = response.read(), response.status, response.headers
                break
            except urllib.error.HTTPError as error:
                if error.code == 404:
                    raw, code, headers = error.read(), 404, error.headers
                    break
                if error.code not in (429, 500, 502, 503, 504) or attempt == 3:
                    raise CatalogueHTTPError(url, error) from error
                retry = error.headers.get('Retry-After', '')
                time.sleep(min(30, int(retry)) if retry.isdigit() else 2 ** attempt)
            except (urllib.error.URLError, TimeoutError):
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)
        I.require(code in (200, 404), 'Unexpected catalogue HTTP status')
        meta = dict(url=url, sha256=I.sha(raw), http_status=code,
                    retrieved_utc=datetime.now(timezone.utc).isoformat(),
                    service_release=headers.get('X-UniProt-Release', 'not_supplied_by_service'))
        # Validate successful JSON before caching; cached pairs retain raw error bodies for 404.
        if code == 200:
            json.loads(raw)
        cache.mkdir(parents=True, exist_ok=True)
        body.write_bytes(raw)
        sidecar.write_text(json.dumps(meta, indent=2) + '\n')
    return (json.loads(raw) if meta['http_status'] == 200 else None), meta


def reference_cache(project, binding, sequences):
    """Import frozen public response bytes into a separate, verified namespace."""
    descriptor = project / 'references/step05a_metadata_2026-10-09/reference.json'
    ref = json.loads(descriptor.read_text())
    I.require(ref['protocol_version'] == '0.1.0' and
              ref['purpose'] == 'frozen_public_API_response_reference_not_workstation_download' and
              ref['canonical_fasta_sha256'] == binding['input_sha256'][binding['fasta_file']],
              'Reference protocol/purpose/sequence universe mismatch')
    archive = I.safe_path(project, ref['archive_file'])
    ledger = I.safe_path(project, ref['source_manifest_file'])
    I.require(I.sha(archive.read_bytes()) == ref['archive_sha256'] and
              I.sha(ledger.read_bytes()) == ref['source_manifest_sha256'], 'Reference checksum mismatch')
    records = I.read_csv(ledger)
    urls = {base + accession for base in (PDBE, AFDB) for accession in sequences}
    I.require(len(records) == ref['query_records'] == len(urls) and
              {r['url'] for r in records} == urls, 'Reference query universe mismatch')
    names = {I.sha(url.encode()) + suffix for url in urls for suffix in ('.response', '.metadata.json')}
    contents = {}
    with tarfile.open(archive, 'r:gz') as bundle:
        for member in bundle:
            I.require(member.isfile() and member.name in names and member.name not in contents,
                      'Unsafe/unexpected/duplicate reference archive member')
            I.require(0 <= member.size <= ref['uncompressed_bytes'], 'Reference member size invalid')
            contents[member.name] = bundle.extractfile(member).read()
    I.require(set(contents) == names and len(contents) == ref['cache_files'] and
              sum(map(len, contents.values())) == ref['uncompressed_bytes'], 'Reference archive incomplete')
    for row in records:
        key = I.sha(row['url'].encode())
        raw = contents[key + '.response']
        meta = json.loads(contents[key + '.metadata.json'])
        I.require(meta['url'] == row['url'] and meta['sha256'] == row['sha256'] == I.sha(raw) and
                  meta['http_status'] == int(row['http_status']) and meta['http_status'] in (200, 404) and
                  meta['retrieved_utc'] == row['retrieved_utc'] and
                  meta['service_release'] == row['service_release'], 'Reference response provenance mismatch')
        if meta['http_status'] == 200:
            json.loads(raw)
    cache = project / 'data/raw/structure_catalogue/reference' / ref['archive_sha256']
    if cache.exists():
        I.require({p.name for p in cache.iterdir()} == names and
                  all((cache / name).read_bytes() == raw for name, raw in contents.items()),
                  'Existing reference cache changed; remove that reference cache directory before retrying')
    else:
        cache.parent.mkdir(parents=True, exist_ok=True)
        stage = Path(tempfile.mkdtemp(prefix='reference_stage_', dir=cache.parent))
        try:
            for name, raw in contents.items():
                (stage / name).write_bytes(raw)
            stage.rename(cache)
        finally:
            if stage.exists():
                shutil.rmtree(stage)
    return cache, ref, [descriptor, archive, ledger]


def number(value, label):
    I.require(not isinstance(value, bool), 'Boolean in numeric field: ' + label)
    n = float(value)
    I.require(math.isfinite(n) and n.is_integer(), 'Invalid integer: ' + label)
    return int(n)


def current_field(record, current, legacy):
    if current in record:
        I.require(legacy not in record or record[legacy] == record[current], 'Conflicting AlphaFold schema fields')
        return record[current]
    I.require(legacy in record, 'Missing AlphaFold field: ' + current)
    return record[legacy]


def parse_pdbe(accession, data, canonical_length):
    if data is None:
        return []
    I.require(isinstance(data, dict) and accession in data and isinstance(data[accession], list),
              'Unexpected PDBe schema for ' + accession)
    records, seen = [], set()
    for item in data[accession]:
        pdb, chain = item['pdb_id'], item['chain_id']
        I.require(isinstance(pdb, str) and re.fullmatch(r'[A-Za-z0-9_]+', pdb) and
                  isinstance(chain, str) and bool(chain.strip()), 'Invalid PDB/chain identifier')
        start, end = number(item['unp_start'], 'unp_start'), number(item['unp_end'], 'unp_end')
        I.require(1 <= start <= end, 'Invalid coarse PDBe UniProt interval')
        taxid = number(item['tax_id'], 'tax_id') if item.get('tax_id') is not None else ''
        signature = json.dumps(item, sort_keys=True)
        if signature in seen:
            continue
        seen.add(signature)
        if taxid == 10090:
            status, reason = 'candidate', 'MOUSE_SIFTS_CANDIDATE_LOCAL_REVIEW_PENDING'
        elif taxid == '':
            status, reason = 'held', 'PDB_CHAIN_TAXONOMY_NOT_PROVIDED'
        else:
            status, reason = 'excluded', 'NON_MOUSE_PDB_CHAIN'
        if taxid == 10090 and end > canonical_length:
            status, reason = 'held', 'PDBE_INTERVAL_OUTSIDE_FROZEN_CANONICAL'
        records.append(dict(protein_accession=accession, source='PDBe_SIFTS', structure_id=pdb,
                            chain_id=chain, model_version='', species_taxid=taxid,
                            sequence_start=start, sequence_end=end, sequence_match='not_assessed',
                            model_sequence_sha256='', reported_sequence_checksum='',
                            experimental_method=item.get('experimental_method', ''),
                            resolution=item.get('resolution', ''), reported_coverage=item.get('coverage', ''),
                            cif_url='', pdb_url='', confidence_url='', pae_url='',
                            candidate_status=status, reason=reason))
    return records


def parse_afdb(accession, canonical, data):
    if data is None:
        return []
    I.require(isinstance(data, list), 'Unexpected AlphaFold schema for ' + accession)
    records, seen = [], set()
    for item in data:
        identifier = current_field(item, 'modelEntityId', 'entryId')
        sequence = current_field(item, 'sequence', 'uniprotSequence')
        start = number(current_field(item, 'sequenceStart', 'uniprotStart'), 'sequenceStart')
        end = number(current_field(item, 'sequenceEnd', 'uniprotEnd'), 'sequenceEnd')
        I.require(isinstance(identifier, str) and isinstance(sequence, str) and
                  bool(sequence) and set(sequence) <= I.ALPHABET and 1 <= start <= end,
                  'Invalid AlphaFold model/sequence metadata')
        version = number(item['latestVersion'], 'latestVersion')
        I.require(version >= 1, 'Invalid AlphaFold version')
        taxid = number(item['taxId'], 'taxId') if item.get('taxId') is not None else ''
        model_accession = item.get('uniprotAccession', '')
        if model_accession != accession:
            status, reason, match = 'held', 'ALPHAFOLD_ACCESSION_DIFFERS_FROM_CANONICAL', 'not_assessed'
        elif taxid != 10090:
            status, reason, match = 'held', 'ALPHAFOLD_MOUSE_TAXONOMY_NOT_CONFIRMED', 'not_assessed'
        elif item.get('isComplex', False):
            status, reason, match = 'held', 'ALPHAFOLD_COMPLEX_REQUIRES_SEPARATE_MAPPING', 'not_assessed'
        elif end > len(canonical):
            status, reason, match = 'held', 'ALPHAFOLD_INTERVAL_OUTSIDE_CANONICAL', 'mismatch'
        elif sequence == canonical and start == 1 and end == len(canonical):
            status, reason, match = 'candidate', 'EXACT_FULL_SEQUENCE_LOCAL_QUALITY_PENDING', 'exact_full'
        elif sequence == canonical[start - 1:end]:
            status, reason, match = 'candidate', 'EXACT_PARTIAL_SEQUENCE_BOUNDARY_REVIEW_PENDING', 'exact_partial'
        else:
            status, reason, match = 'held', 'ALPHAFOLD_SEQUENCE_DIFFERS_FROM_FROZEN_CANONICAL', 'mismatch'
        urls = {out: item.get(key, '') for out, key in
                [('cif_url', 'cifUrl'), ('pdb_url', 'pdbUrl'), ('confidence_url', 'plddtDocUrl'), ('pae_url', 'paeDocUrl')]}
        for value in urls.values():
            if value:
                parts = urllib.parse.urlsplit(value)
                I.require(parts.scheme == 'https' and parts.netloc in ('alphafold.ebi.ac.uk', 'alphafold.com',
                                                                      'www.alphafold.ebi.ac.uk'), 'Unexpected AlphaFold file URL')
        if status == 'candidate' and not urls['cif_url'] and not urls['pdb_url']:
            status, reason = 'held', 'ALPHAFOLD_COORDINATE_URL_NOT_PROVIDED'
        signature = json.dumps(item, sort_keys=True)
        if signature in seen:
            continue
        seen.add(signature)
        records.append(dict(protein_accession=accession, source='AlphaFold_DB', structure_id=identifier,
                            chain_id=item.get('chainId', ''), model_version=version, species_taxid=taxid,
                            sequence_start=start, sequence_end=end, sequence_match=match,
                            model_sequence_sha256=I.sha(sequence.encode()),
                            reported_sequence_checksum=item.get('sequenceChecksum', ''),
                            experimental_method='predicted', resolution='', reported_coverage='',
                            candidate_status=status, reason=reason, **urls))
    return records


def catalogue_one(accession, sequence, cache, offline):
    pdbe, p_source = fetch(PDBE + accession, cache, offline)
    afdb, a_source = fetch(AFDB + accession, cache, offline)
    candidates = parse_pdbe(accession, pdbe, len(sequence)) + parse_afdb(accession, sequence, afdb)
    for row in candidates:
        source = p_source if row['source'] == 'PDBe_SIFTS' else a_source
        row.update(query_url=source['url'], response_sha256=source['sha256'],
                   structural_eligibility='not_assessed', protocol_version='0.1.0')
    p_count = sum(r['source'] == 'PDBe_SIFTS' and r['candidate_status'] == 'candidate' for r in candidates)
    a_count = sum(r['source'] == 'AlphaFold_DB' and r['candidate_status'] == 'candidate' for r in candidates)
    if p_count or a_count:
        status, reason = 'held', 'STRUCTURAL_CANDIDATES_REQUIRE_LOCAL_REVIEW'
    elif any(r['candidate_status'] == 'held' for r in candidates):
        status, reason = 'held', 'ONLY_UNRESOLVED_METADATA_CANDIDATES_REPORTED'
    else:
        status, reason = 'unassessable', 'NO_METADATA_ELIGIBLE_CANDIDATE_REPORTED'
    inventory = dict(protein_accession=accession, sequence_length=len(sequence), sequence_checksum=I.sha(sequence.encode()),
                     pdbe_records=sum(r['source'] == 'PDBe_SIFTS' for r in candidates), mouse_pdbe_candidates=p_count,
                     alphafold_records=sum(r['source'] == 'AlphaFold_DB' for r in candidates),
                     exact_sequence_alphafold_candidates=a_count, status=status,
                     reason=reason, protocol_version='0.1.0')
    return candidates, inventory, [dict(p_source, protein_accession=accession, service='PDBe_SIFTS'),
                                   dict(a_source, protein_accession=accession, service='AlphaFold_DB')]


def run(project, offline=False, workers=8, use_reference_cache=False):
    out = project / 'results/05a_structure_catalogue'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'SUCCESS.txt').unlink(missing_ok=True)
    (out / 'failure_context.json').unlink(missing_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='catalogue_stage_', dir=out))
    try:
        I.require(1 <= workers <= 16, 'Workers must be between 1 and 16')
        binding, sequences, inputs = load_input(project)
        cache = project / 'data/raw/structure_catalogue/step05a' / binding['input_sha256'][binding['fasta_file']]
        reference = None
        if use_reference_cache:
            cache, reference, reference_inputs = reference_cache(project, binding, sequences)
            inputs += reference_inputs
            offline = True
        cache.mkdir(parents=True, exist_ok=True)
        candidates, proteins, sources = [], [], []
        print(f'[1/3] Cataloguing PDBe and AlphaFold metadata for {len(sequences)} proteins...', flush=True)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            jobs = [pool.submit(catalogue_one, a, seq, cache, offline) for a, seq in sorted(sequences.items())]
            try:
                for n, job in enumerate(as_completed(jobs), 1):
                    records, protein, provenance = job.result()
                    candidates.extend(records); proteins.append(protein); sources.extend(provenance)
                    if n % 25 == 0 or n == len(jobs):
                        print(f'  Proteins complete: {n}/{len(jobs)}', flush=True)
            except Exception:
                for job in jobs:
                    job.cancel()
                raise
        candidates.sort(key=lambda r: (r['protein_accession'], r['source'], r['structure_id'],
                                       str(r['chain_id']), r['sequence_start'], r['sequence_end']))
        proteins.sort(key=lambda r: r['protein_accession']); sources.sort(key=lambda r: r['url'])
        print('[2/3] Saving candidates and complete protein/query audits; no structure selected...', flush=True)
        candidate_fields = ['protein_accession', 'source', 'structure_id', 'chain_id', 'model_version', 'species_taxid',
                            'sequence_start', 'sequence_end', 'sequence_match', 'model_sequence_sha256',
                            'reported_sequence_checksum', 'experimental_method', 'resolution', 'reported_coverage',
                            'cif_url', 'pdb_url', 'confidence_url', 'pae_url', 'candidate_status', 'reason',
                            'query_url', 'response_sha256', 'structural_eligibility', 'protocol_version']
        I.write_csv(stage / 'structure_candidates.csv', candidates, candidate_fields)
        I.write_csv(stage / 'protein_structure_inventory.csv', proteins, list(proteins[0]))
        I.write_csv(stage / 'source_manifest.csv', sources, list(sources[0]))
        audit = [dict(step_id='05', object_level='protein', object_id=r['protein_accession'], status=r['status'],
                      reason_code=r['reason'], reason_text=r['reason'], input_source=binding['fasta_file'],
                      input_checksum=binding['input_sha256'][binding['fasta_file']], protocol_version='0.1.0') for r in proteins]
        I.write_csv(stage / 'step_audit.csv', audit, list(audit[0]))
        inputs = list(dict.fromkeys(inputs + [Path(__file__).resolve(), Path(I.__file__).resolve()]))
        I.write_csv(stage / 'input_checksums.csv', [dict(input=str(p.relative_to(project)), sha256=I.sha(p.read_bytes()))
                                                   for p in inputs], ['input', 'sha256'])
        summary = dict(status='STEP05A_GENERATED_REVIEW_PENDING', implementation_version=VERSION, protocol_version='0.1.0',
                       input_proteins=len(proteins), input_sites=binding['input_sites'], candidate_records=len(candidates),
                       query_records=len(sources), proteins_with_mouse_pdbe_candidates=sum(r['mouse_pdbe_candidates'] > 0 for r in proteins),
                       proteins_with_exact_sequence_alphafold=sum(r['exact_sequence_alphafold_candidates'] > 0 for r in proteins),
                       proteins_with_both=sum(r['mouse_pdbe_candidates'] > 0 and r['exact_sequence_alphafold_candidates'] > 0 for r in proteins),
                       proteins_without_metadata_eligible_candidate=sum(not r['mouse_pdbe_candidates'] and
                                                                         not r['exact_sequence_alphafold_candidates'] for r in proteins),
                       protein_status_counts=dict(Counter(r['status'] for r in proteins)),
                       candidate_status_counts=dict(Counter(r['candidate_status'] for r in candidates)),
                       candidate_reason_counts=dict(Counter(r['reason'] for r in candidates)),
                       alphafold_model_versions=sorted({r['model_version'] for r in candidates if r['source'] == 'AlphaFold_DB'}),
                       http_status_counts=dict(Counter(str(r['http_status']) for r in sources)),
                       metadata_mode='frozen_reference' if reference else 'workstation_cache_or_network',
                       reference_bundle_sha256=reference['archive_sha256'] if reference else '',
                       workers=workers, offline=offline, python=platform.python_version(), platform=platform.platform(),
                       generated_utc=datetime.now(timezone.utc).isoformat())
        (stage / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
        report = ['# Step 05a — structural candidate catalogue', '', 'Status: ' + summary['status'], '',
                  f'- Input proteins: {len(proteins)}; canonical Cys sites: {binding["input_sites"]}',
                  f'- Proteins with mouse PDBe candidates: {summary["proteins_with_mouse_pdbe_candidates"]}',
                  f'- Proteins with sequence-matching AlphaFold candidates: {summary["proteins_with_exact_sequence_alphafold"]}',
                  f'- Proteins with both sources: {summary["proteins_with_both"]}',
                  f'- Proteins without metadata-eligible candidates: {summary["proteins_without_metadata_eligible_candidate"]}',
                  f'- Candidate records: {len(candidates)}; completed queries: {len(sources)}', '',
                  'This is a candidate catalogue within Step 05, not structural inclusion approval.',
                  'Every offered PDBe chain and AlphaFold record is audited; API order does not choose a structure.',
                  'PDBe coarse intervals are not residue alignments or proof of observed Cys/SG atoms.',
                  'AlphaFold identity checks use frozen canonical sequence, taxid and accession; partial models retain boundary-review flags.',
                  'No coordinate files, assemblies, mutations, atom completeness or local quality have been assessed.',
                  'Target pLDDT >=90 and neighbor pLDDT >=70 within 6 A remain fixed for the next substep.',
                  '404/empty query results mean not reported by that service, not absent protein or oxidative resistance.',
                  'Network/server errors fail the run and are never counted as missing structures.',
                  'Raw response hashes and retrieval times are cached separately from model versions; no global current database release is invented.',
                  'Review the catalogue before coordinate retrieval, per-site mapping and quality assessment in Step 05b.']
        if reference:
            report += ['', 'Metadata source: frozen public API reference captured on 2026-10-09; this run made no API requests.',
                       'Original upstream retrieval times and response hashes are preserved; these are not fresh workstation downloads.',
                       'Reference bundle SHA256: ' + reference['archive_sha256']]
        (stage / 'catalogue_report.md').write_text('\n'.join(report) + '\n')
        for path in stage.iterdir():
            shutil.copyfile(path, out / path.name)
        (out / 'FAILURE.txt').unlink(missing_ok=True)
        (out / 'SUCCESS.txt').write_text(summary['status'] + '\n')
        print('[3/3] Outputs: ' + str(out), flush=True)
        print('\n'.join(report[:10]))
        return summary
    except Exception as error:
        (out / 'FAILURE.txt').write_text(str(error) + '\n')
        if isinstance(error, CatalogueHTTPError):
            (out / 'failure_context.json').write_text(json.dumps(error.details, indent=2) + '\n')
        raise
    finally:
        shutil.rmtree(stage)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-dir', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--offline', action='store_true', help='Require complete checksum-verified metadata caches')
    parser.add_argument('--reference-cache', action='store_true',
                        help='Verify and use the bundled frozen API reference without network requests')
    parser.add_argument('--workers', type=int, default=8, help='Concurrent protein queries, 1–16 (default 8)')
    args = parser.parse_args()
    try:
        run(args.project_dir.resolve(), args.offline, args.workers, args.reference_cache)
    except Exception as error:
        print('STEP05A FAILED: ' + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
