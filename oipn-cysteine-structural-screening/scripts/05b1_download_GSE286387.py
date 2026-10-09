#!/usr/bin/env python3
"""Step 05b1: resumable raw coordinate acquisition/archive; no structural approval."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request

SPEC = importlib.util.spec_from_file_location('catalogue', Path(__file__).with_name('05a_catalogue_GSE286387.py'))
C = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(C)
I = C.I
VERSION = '1.0.0'
HOSTS = {'files.rcsb.org', 'alphafold.ebi.ac.uk'}
KINDS = ('pdb_mmcif', 'alphafold_mmcif', 'alphafold_confidence', 'alphafold_pae')


def utc():
    return datetime.now(timezone.utc).isoformat()


def permitted(url):
    u = urllib.parse.urlsplit(url)
    I.require(u.scheme == 'https' and u.hostname in HOSTS and not u.username and not u.password
              and u.port in (None, 443), 'Unapproved structure download URL: ' + url)


def load_input(project):
    path = project / 'config/GSE286387_step05b_input.json'
    binding = json.loads(path.read_text())
    I.require(binding['status'] == 'STEP05A_SNAPSHOT_ACCEPTED_FOR_STEP05B' and
              binding['species_taxid'] == 10090 and binding['protocol_version'] == '0.1.0',
              'Unaccepted Step 05b input')
    inputs = [path]
    for relative, checksum in binding['input_sha256'].items():
        f = I.safe_path(project, relative)
        I.require(f.is_file() and I.sha(f.read_bytes()) == checksum, 'Frozen Step 05b input changed: ' + relative)
        inputs.append(f)
    for key in ('candidate_file', 'protein_file', 'source_file'):
        I.require(binding[key] in binding['input_sha256'], 'Unbound primary catalogue input')
    previous, sequences, upstream = C.load_input(project)
    I.require(binding['input_proteins'] == len(sequences) == 760 and
              binding['input_sites'] == previous['input_sites'] == 10799, 'Canonical input universe mismatch')
    snapshot = I.read_csv(I.safe_path(project, binding['candidate_file']))
    recorded = I.read_csv(I.safe_path(project, binding['candidate_file']).parent / 'input_checksums.csv')
    for row in recorded:
        f = I.safe_path(project, row['input'])
        I.require(I.sha(f.read_bytes()) == row['sha256'], 'Catalogue provenance input changed: ' + row['input'])
        inputs.append(f)
    I.require(all(r['protein_accession'] in sequences and r['structural_eligibility'] == 'not_assessed'
                  and r['protocol_version'] == '0.1.0' for r in snapshot), 'Catalogue identity/eligibility changed')
    return binding, snapshot, list(dict.fromkeys(inputs + upstream + [Path(__file__).resolve()]))


def make_plan(rows):
    files, links = {}, []
    for n, row in enumerate(rows, 1):
        ids = []
        if row['candidate_status'] == 'candidate':
            I.require(row['species_taxid'] == '10090', 'Nonmouse metadata candidate')
            if row['source'] == 'PDBe_SIFTS':
                name = row['structure_id'].lower()
                I.require(re.fullmatch(r'[a-z0-9_]+', name) is not None, 'Unsafe PDB identifier')
                options = [('pdb_mmcif', 'https://files.rcsb.org/download/' + name + '.cif', name, True)]
            elif row['source'] == 'AlphaFold_DB':
                I.require(row['sequence_match'] in ('exact_full', 'exact_partial'), 'Unmatched predicted candidate')
                options = [('alphafold_mmcif', row['cif_url'], row['structure_id'], True)]
                options += [(kind, row[field], row['structure_id'], False) for kind, field in
                            [('alphafold_confidence', 'confidence_url'), ('alphafold_pae', 'pae_url')] if row[field]]
            else:
                raise ValueError('Unknown catalogue source')
            for kind, url, structure, required in options:
                I.require(bool(url), 'Required coordinate URL missing')
                permitted(url)
                fid = I.sha((kind + '\n' + url).encode())
                item = dict(file_id=fid, kind=kind, url=url, structure_id=structure, required=required)
                if fid in files:
                    I.require(files[fid] == item, 'Conflicting duplicate file metadata')
                files[fid] = item; ids.append(fid)
        links.append(dict(candidate_row=n, protein_accession=row['protein_accession'], source=row['source'],
                          structure_id=row['structure_id'], chain_id=row['chain_id'],
                          candidate_status=row['candidate_status'], reason=row['reason'],
                          file_ids=';'.join(ids), acquisition_action='planned' if ids else 'metadata_held_or_excluded',
                          structural_eligibility='not_assessed', protocol_version='0.1.0'))
    return sorted(files.values(), key=lambda x: (KINDS.index(x['kind']), x['url'])), links


def validate_payload(item, raw):
    I.require(bool(raw), 'Empty downloaded file')
    if item['kind'].endswith('mmcif'):
        text = raw.decode('utf-8')
        block = re.search(r'^data_(\S+)', text, re.M)
        I.require(block is not None and block.group(1).casefold() == item['structure_id'].casefold(),
                  'mmCIF data-block identifier mismatch')
        I.require(all(re.search(r'^' + re.escape(tag) + r'(?:\s|$)', text, re.M) for tag in
                      ('_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z')),
                  'Downloaded file lacks coordinate tags')
    else:
        value = json.loads(raw)
        I.require(isinstance(value, (dict, list)) and bool(value), 'Empty/noncontainer auxiliary JSON')
    # Transport/content sanity only. Full parsing, sequence/atom mapping and quality are future work.


def paths(cache, item):
    suffix = '.cif' if item['kind'].endswith('mmcif') else '.json'
    return cache / (item['file_id'] + suffix), cache / (item['file_id'] + '.metadata.json')


def cached(cache, item):
    body, sidecar = paths(cache, item)
    if not sidecar.exists():
        return None  # Uncommitted/orphan bytes are never accepted; a fresh download will replace them.
    I.require(body.is_file(), 'Cached structure body missing: ' + item['url'])
    raw, meta = body.read_bytes(), json.loads(sidecar.read_text())
    I.require(meta['url'] == item['url'] and meta['file_id'] == item['file_id'] and
              meta['sha256'] == I.sha(raw) and meta['bytes'] == len(raw) and meta['http_status'] == 200,
              'Corrupt structure cache: ' + item['url'])
    permitted(meta['response_url']); validate_payload(item, raw)
    return {**item, **meta, 'status': 'downloaded', 'reason': 'VERIFIED_CACHE_REUSED',
            'local_file': str(body), 'metadata_file': str(sidecar)}


def atomic(path, raw):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
        temporary = Path(stream.name); stream.write(raw)
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def download(cache, item, offline=False, timeout=30, attempts=3):
    try:
        result = cached(cache, item)
        if result:
            return result
        I.require(not offline, 'Offline structure cache missing')
        for attempt in range(attempts):
            try:
                request = urllib.request.Request(item['url'], headers={'User-Agent': 'OIPN-coordinate-archive/' + VERSION,
                                                                      'Accept-Encoding': 'identity'})
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    permitted(response.geturl())
                    I.require(response.status == 200, 'Unexpected coordinate HTTP status')
                    raw = response.read()
                    declared = response.headers.get('Content-Length')
                    if declared is not None:
                        I.require(len(raw) == int(declared), 'Incomplete downloaded payload')
                    meta = dict(file_id=item['file_id'], url=item['url'], response_url=response.geturl(),
                                http_status=200, bytes=len(raw), sha256=I.sha(raw), retrieved_utc=utc(),
                                etag=response.headers.get('ETag', ''), last_modified=response.headers.get('Last-Modified', ''))
                validate_payload(item, raw)
                body, sidecar = paths(cache, item)
                atomic(body, raw)
                atomic(sidecar, (json.dumps(meta, indent=2) + '\n').encode())
                return {**item, **meta, 'status': 'downloaded', 'reason': 'DOWNLOADED_CONTENT_CHECKED',
                        'local_file': str(body), 'metadata_file': str(sidecar)}
            except urllib.error.HTTPError as error:
                if error.code == 404:
                    return dict(item, status='unavailable', reason='HTTP_404_NOT_REPORTED', http_status=404,
                                retrieved_utc=utc(), error=str(error)[:1000])
                if error.code not in (429, 500, 502, 503, 504) or attempt == attempts - 1:
                    raise
                delay = error.headers.get('Retry-After', '') if error.headers else ''
                time.sleep(min(30, int(delay)) if delay.isdigit() else 2 ** attempt)
            except (urllib.error.URLError, TimeoutError, ConnectionError, OSError):
                if attempt == attempts - 1:
                    raise
                time.sleep(2 ** attempt)
    except Exception as error:
        return dict(item, status='failed', reason=type(error).__name__, retrieved_utc=utc(), error=str(error)[:1000])


def archive_files(project, out, cache, items, inputs, plan_hash):
    files = list(dict.fromkeys(inputs + [f for f in out.iterdir() if f.is_file() and f.name not in ('archive_receipt.json', 'archive_manifest.csv')]))
    for item in items:
        if item['status'] == 'downloaded':
            I.require(cached(cache, item) is not None, 'Archive cache verification failed')
            files += list(paths(cache, item))
    files = sorted(set(files))
    manifest = [dict(path=str(f.relative_to(project)), sha256=I.sha(f.read_bytes()), bytes=f.stat().st_size) for f in files]
    I.write_csv(out / 'archive_manifest.csv', manifest, ['path', 'sha256', 'bytes'])
    files.append(out / 'archive_manifest.csv')
    destination = project / 'data/archives' / ('step05b1_' + plan_hash[:16] + '.tar.gz')
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix('.partial')
    try:
        with tarfile.open(temporary, 'w:gz') as bundle:
            for f in files:
                I.require(f.is_file() and not f.is_symlink(), 'Nonregular archive input')
                bundle.add(f, arcname=str(f.relative_to(project)), recursive=False)
        expected = {str(f.relative_to(project)): I.sha(f.read_bytes()) for f in files}
        with tarfile.open(temporary, 'r:gz') as bundle:
            found = {}
            for member in bundle:
                I.require(member.isfile() and member.name in expected and member.name not in found, 'Archive member mismatch')
                digest = hashlib.sha256()
                with bundle.extractfile(member) as stream:
                    while block := stream.read(1024 * 1024):
                        digest.update(block)
                found[member.name] = digest.hexdigest()
            I.require(found == expected, 'Archive content checksum verification failed')
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    digest = hashlib.sha256()
    with destination.open('rb') as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    receipt = dict(archive=str(destination.relative_to(project)), sha256=digest.hexdigest(),
                   bytes=destination.stat().st_size, members=len(files), generated_utc=utc(),
                   scope='downloaded_raw_files_and_provenance_not_structural_approval')
    atomic(out / 'archive_receipt.json', (json.dumps(receipt, indent=2) + '\n').encode())
    return receipt


def run(project, workers=4, offline=False, probe=False, plan_only=False, archive=False):
    out = project / 'results/05b1_structure_download'
    out.mkdir(parents=True, exist_ok=True)
    for marker in ('SUCCESS.txt', 'PROBE_SUCCESS.txt', 'FAILURE.txt'):
        (out / marker).unlink(missing_ok=True)
    (out / 'archive_receipt.json').unlink(missing_ok=True)
    try:
        I.require(1 <= workers <= 8, 'Workers must be between 1 and 8')
        I.require(not (archive and (probe or plan_only)), '--archive requires a full acquisition run')
        binding, rows, inputs = load_input(project)
        plan, links = make_plan(rows)
        plan_hash = I.sha(json.dumps(plan, sort_keys=True).encode())
        cache = project / 'data/raw/structure_coordinates/step05b1' / plan_hash
        I.write_csv(out / 'download_plan.csv', plan, list(plan[0]))
        I.write_csv(out / 'candidate_file_links.csv', links, list(links[0]))
        I.write_csv(out / 'input_checksums.csv', [dict(input=str(f.relative_to(project)), sha256=I.sha(f.read_bytes()))
                                                 for f in inputs], ['input', 'sha256'])
        print('File plan: ' + str(dict(Counter(r['kind'] for r in plan))), flush=True)
        if plan_only:
            print('PLAN_ONLY: no files downloaded; no acquisition approval.', flush=True)
            return 0
        # Check one endpoint per file kind before any bulk requests, including on resumption.
        first = {kind: next(r for r in plan if r['kind'] == kind) for kind in KINDS if any(r['kind'] == kind for r in plan)}
        results = {r['file_id']: download(cache, r, offline) for r in first.values()}
        blocked = any(r['status'] == 'failed' or (r['required'] and r['status'] != 'downloaded') for r in results.values())
        if not probe and not blocked:
            pending = [r for r in plan if r['file_id'] not in results]
            with ThreadPoolExecutor(max_workers=workers) as pool:
                jobs = {pool.submit(download, cache, r, offline): r for r in pending}
                for n, job in enumerate(as_completed(jobs), len(results) + 1):
                    r = job.result(); results[r['file_id']] = r
                    if n % 25 == 0 or n == len(plan):
                        print(f'  Files checked: {n}/{len(plan)}', flush=True)
        records = [results.get(r['file_id'], dict(r, status='pending', reason='PROBE_ONLY' if probe else 'PREFLIGHT_BLOCKED')) for r in plan]
        fields = list(plan[0]) + ['status', 'reason', 'http_status', 'bytes', 'sha256', 'retrieved_utc',
                                 'response_url', 'etag', 'last_modified', 'local_file', 'metadata_file', 'error']
        for r in records:
            for key in ('local_file', 'metadata_file'):
                if r.get(key): r[key] = str(Path(r[key]).relative_to(project))
        I.write_csv(out / 'download_manifest.csv', records, fields)
        required = [r for r in records if r['required']]
        complete = all(r['status'] == 'downloaded' for r in required)
        errors = sum(r['status'] == 'failed' for r in records)
        status = ('PROBE_COMPLETE' if not blocked else 'PROBE_FAILED') if probe else (
            'STEP05B1_GENERATED_REVIEW_PENDING' if complete and not errors else 'STEP05B1_INCOMPLETE')
        summary = dict(status=status, implementation_version=VERSION, protocol_version='0.1.0',
                       source_commit=binding['source_commit'], input_proteins=760, input_sites=10799,
                       offered_candidates=len(rows), file_plan_sha256=plan_hash, files_planned=len(plan),
                       file_kind_counts=dict(Counter(r['kind'] for r in plan)),
                       file_status_counts=dict(Counter(r['status'] for r in records)),
                       required_coordinates=len(required), coordinates_downloaded=sum(r['status']=='downloaded' for r in required),
                       auxiliary_unavailable=sum(not r['required'] and r['status']=='unavailable' for r in records),
                       offline=offline, probe=probe, workers=workers, generated_utc=utc())
        atomic(out / 'summary.json', (json.dumps(summary, indent=2) + '\n').encode())
        report = '\n'.join(['# Step 05b1 — raw structure acquisition', '', 'Status: ' + status, '',
                           f'- Required coordinates: {summary["coordinates_downloaded"]}/{len(required)}',
                           '- File statuses: ' + str(summary['file_status_counts']),
                           '- Auxiliary files unavailable (404): ' + str(summary['auxiliary_unavailable']), '',
                           'Raw mmCIF is preserved; PDB entries are complete deposited asymmetric units, not selected chains or generated biological assemblies.',
                           'No residue mapping, native-state/assembly/sequence validation, local validation, pLDDT gate or structural pass has been issued.',
                           'Only metadata candidates are downloaded; all held/excluded offered records remain in candidate_file_links.csv.',
                           'Network errors are acquisition failures, never biological absence. Auxiliary 404s remain explicit evidence gaps.',
                           'Caches reuse only matching byte hashes and original retrieval provenance. No TLS verification is disabled.',
                           'Archive completion refers to raw coordinate acquisition only; inspect manifests before later structural work.']) + '\n'
        atomic(out / 'download_report.md', report.encode())
        failures = [r for r in records if r['status'] in ('failed', 'unavailable')]
        atomic(out / 'failure_context.json', (json.dumps(failures, indent=2) + '\n').encode())
        ok = not blocked if probe else complete and not errors
        marker = 'PROBE_SUCCESS.txt' if probe and ok else 'SUCCESS.txt' if ok else 'FAILURE.txt'
        atomic(out / marker, (status + '\n').encode())
        if archive:
            # Pending/error cases remain in the archive manifests; the archive does not hide partial acquisition.
            archive_files(project, out, cache, records, inputs, plan_hash)
        print(report, flush=True)
        print('Outputs: ' + str(out), flush=True)
        return 0 if ok else 1
    except Exception as error:
        (out / 'SUCCESS.txt').unlink(missing_ok=True)
        (out / 'PROBE_SUCCESS.txt').unlink(missing_ok=True)
        atomic(out / 'FAILURE.txt', (str(error) + '\n').encode())
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-dir', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--offline', action='store_true')
    parser.add_argument('--probe', action='store_true', help='Test one file per kind; no full-download success marker')
    parser.add_argument('--plan-only', action='store_true')
    parser.add_argument('--archive', action='store_true', help='Create and verify a tar.gz of available raw files and provenance')
    args = parser.parse_args()
    try:
        return run(args.project_dir.resolve(), args.workers, args.offline, args.probe, args.plan_only, args.archive)
    except Exception as error:
        print('STEP05B1 FAILED: ' + str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
