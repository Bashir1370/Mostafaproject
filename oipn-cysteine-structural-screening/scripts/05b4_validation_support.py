"""Resumable, provenance-preserving official wwPDB XML acquisition; no quality cutoffs."""
from datetime import datetime, timezone
import gzip
import hashlib
import http.client
import json
from pathlib import Path
import re
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

VERSION = '1.0.0'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def utc():
    return datetime.now(timezone.utc).isoformat()


def permitted(url):
    u = urllib.parse.urlsplit(url)
    require(u.scheme == 'https' and u.hostname == 'files.wwpdb.org' and
            u.path.startswith('/pub/pdb/validation_reports/'), 'Unexpected validation URL: ' + url)


def plan(ids):
    out = []
    for pdb in sorted(set(ids)):
        require(re.fullmatch(r'[0-9][a-z0-9]{3}', pdb) is not None, 'Unsupported frozen PDB identifier')
        url = f'https://files.wwpdb.org/pub/pdb/validation_reports/{pdb[1:3]}/{pdb}/{pdb}_validation.xml.gz'
        out.append(dict(file_id=sha(url.encode()), structure_id=pdb, url=url))
    return out


def xml_content(raw, pdb):
    root = ET.fromstring(gzip.decompress(raw))
    entries = [x for x in root.iter() if x.tag.split('}')[-1] == 'Entry']
    require(len(entries) == 1 and entries[0].get('pdbid', '').casefold() == pdb.casefold(),
            'Validation report entry identity mismatch')
    return root, dict(entries[0].attrib)


def atomic(path, raw):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        temp = Path(f.name)
        f.write(raw)
    try:
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def download(cache, item, network=False, timeout=120, attempts=4):
    base = dict(item, status='not_downloaded', reason='OFFLINE_REPORT_NOT_CACHED', http_status='', bytes='',
                sha256='', retrieved_utc='', response_url='', etag='', last_modified='',
                local_file='', metadata_file='', error_excerpt='')
    body, sidecar = cache / (item['file_id'] + '.xml.gz'), cache / (item['file_id'] + '.metadata.json')
    try:
        if sidecar.is_file():
            require(body.is_file(), 'Validation cached body missing')
            raw, meta = body.read_bytes(), json.loads(sidecar.read_text())
            require(all(meta[k] == item[k] for k in ('file_id', 'structure_id', 'url')) and
                    sha(raw) == meta['sha256'] and len(raw) == meta['bytes'] and meta['http_status'] == 200,
                    'Validation cached hash/provenance mismatch')
            permitted(meta['response_url'])
            xml_content(raw, item['structure_id'])
            return dict(base, **meta, status='downloaded', reason='VERIFIED_CACHE_REUSED',
                        local_file=str(body), metadata_file=str(sidecar))
        if not network:
            return base
        permitted(item['url'])
        for attempt in range(attempts):
            try:
                request = urllib.request.Request(item['url'], headers={'User-Agent': 'OIPN-validation/' + VERSION,
                                                                       'Accept-Encoding': 'identity'})
                with urllib.request.urlopen(request, timeout=timeout) as response:
                    permitted(response.geturl())
                    require(response.status == 200, 'Unexpected validation HTTP status')
                    raw = response.read()
                    length = response.headers.get('Content-Length')
                    require(length is None or len(raw) == int(length), 'Incomplete validation HTTP payload')
                    meta = dict(item, response_url=response.geturl(), http_status=200, bytes=len(raw), sha256=sha(raw),
                                retrieved_utc=utc(), etag=response.headers.get('ETag', ''),
                                last_modified=response.headers.get('Last-Modified', ''))
                xml_content(raw, item['structure_id'])
                atomic(body, raw)
                atomic(sidecar, (json.dumps(meta, indent=2) + '\n').encode())
                return dict(base, **meta, status='downloaded', reason='DOWNLOADED_XML_IDENTITY_CHECKED',
                            local_file=str(body), metadata_file=str(sidecar))
            except urllib.error.HTTPError as e:
                transient = e.code == 429 or 500 <= e.code < 600
                if transient and attempt + 1 < attempts:
                    time.sleep(min(2 ** attempt, 8))
                    continue
                return dict(base, status='unavailable' if e.code == 404 else 'failed', reason='HTTP_' + str(e.code),
                            http_status=e.code, retrieved_utc=utc(), response_url=e.geturl(),
                            error_excerpt=e.read(2048).decode('utf-8', errors='replace'))
            except (urllib.error.URLError, TimeoutError, ConnectionError, OSError,
                    http.client.IncompleteRead, http.client.RemoteDisconnected) as e:
                if attempt + 1 == attempts:
                    raise
                time.sleep(min(2 ** attempt, 8))
    except Exception as e:
        return dict(base, status='failed', reason=type(e).__name__ + ': ' + str(e))
    raise RuntimeError('Validation retry loop exhausted unexpectedly')


def residue_index(root):
    """Author-chain/residue/model/insert IDs, with blank altloc normalized explicitly."""
    out = {}
    for element in root.iter():
        if element.tag.split('}')[-1] != 'ModelledSubgroup':
            continue
        a = dict(element.attrib)
        key = tuple(a.get(k, '').strip() for k in ('model', 'chain', 'resnum', 'icode', 'resname'))
        payload = dict(attributes=a, reported_children=[dict(tag=x.tag.split('}')[-1], attributes=dict(x.attrib))
                                                       for x in element])
        out.setdefault(key, []).append(payload)
    return out


def match_residue(index, model, chain, author_number, insertion, component, alternate):
    options = index.get(tuple(str(v).strip() for v in (model, chain, author_number, insertion, component)), [])
    compatible = [p for p in options if p['attributes'].get('altcode', '').strip() in ('', alternate)]
    if len(compatible) != 1:
        return 'not_reported' if not compatible else 'ambiguous_report_records', None
    return 'matched_for_manual_review', compatible[0]
