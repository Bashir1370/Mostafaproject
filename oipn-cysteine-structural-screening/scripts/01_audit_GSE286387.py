#!/usr/bin/env python3
"""Download and audit official GEO inputs. Python standard library only; no DE fit."""
import argparse
import collections
import csv
import datetime
import decimal
import gzip
import hashlib
import json
from pathlib import Path
import platform
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
BASE = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286387/'
FILES = {'soft': ('soft/GSE286387_family.soft.gz'),
         'counts': ('suppl/GSE286387_DRG_raw_counts.txt.gz')}
ANNOTATIONS = ['ensembl_gene_id', 'entrezgene', 'gene_name',
               'gene_biotype', 'external_gene_source']


def table(path, rows, fields, delimiter=','):
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=delimiter)
        w.writeheader()
        w.writerows(rows)


def parse_soft(path):
    samples, series, current = {}, collections.defaultdict(list), None
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        for line in f:
            if line.startswith('^SAMPLE = '):
                key = line.strip().split(' = ', 1)[1]
                if key in samples:
                    raise ValueError('Repeated GEO sample: ' + key)
                current = samples[key] = collections.defaultdict(list)
            elif line.startswith('^'):
                current = None
            elif line.startswith('!') and ' = ' in line:
                key, value = line.rstrip('\n').split(' = ', 1)
                if key.startswith('!Series_'):
                    series[key].append(value)
                elif current is not None:
                    current[key].append(value)
    if series['!Series_geo_accession'] != ['GSE286387']:
        raise ValueError('Unexpected accession')
    if set(series['!Series_sample_id']) != set(samples) or len(samples) != 10:
        raise ValueError('Series/sample membership mismatch')
    return series, samples


def integer(value):
    x = decimal.Decimal(value)
    if not x.is_finite() or x < 0 or x != x.to_integral_value():
        raise ValueError('Counts must be finite nonnegative integers: ' + value)
    return int(x)


def read_counts(path, titles):
    genes, source_rows, duplicates = {}, [], []
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        header = reader.fieldnames
        if not header or len(header) != len(set(header)):
            raise ValueError('Missing or duplicated column labels')
        columns = [x for x in header if x not in ANNOTATIONS]
        if not set(ANNOTATIONS).issubset(header) or set(columns) != set(titles):
            raise ValueError('Counts columns do not exactly match GEO titles')
        for number, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()):
                raise ValueError('Malformed count row ' + str(number))
            key = row['ensembl_gene_id']
            if not key.startswith('ENSMUSG'):
                raise ValueError('Unexpected mouse gene identifier: ' + key)
            values = tuple(integer(row[c]) for c in columns)
            source_rows.append(dict(source_row=number, **{c: row[c] for c in ANNOTATIONS}))
            if key in genes:
                previous = genes[key]
                if previous['values'] != values or previous['biotype'] != row['gene_biotype']:
                    raise ValueError('Conflicting duplicate gene; do not sum: ' + key)
                duplicates.append({'stable_gene_id': key, 'retained_source_row': previous['row'],
                                   'duplicate_source_row': number, 'action': 'keep_one_identical_count_vector'})
                previous['symbols'].add(row['gene_name'])
                previous['entrez'].add(row['entrezgene'])
            else:
                genes[key] = dict(values=values, row=number, biotype=row['gene_biotype'],
                                  symbols={row['gene_name']}, entrez={row['entrezgene']})
    if not genes:
        raise ValueError('Empty matrix')
    return columns, genes, source_rows, duplicates


def run(raw, out, offline=False):
    raw.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    provenance = []
    paths = {}
    for role, suffix in FILES.items():
        path = raw / Path(suffix).name
        url = BASE + suffix
        if not offline:
            temp = path.with_suffix(path.suffix + '.part')
            try:
                with urllib.request.urlopen(url, timeout=60) as response, temp.open('wb') as f:
                    while chunk := response.read(1024 * 1024):
                        f.write(chunk)
                temp.replace(path)
            finally:
                temp.unlink(missing_ok=True)
        if not path.is_file():
            raise ValueError('Missing input: ' + str(path))
        paths[role] = path
        provenance.append(dict(role=role, source_url=url, local_file=path.name,
                               sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                               bytes=path.stat().st_size,
                               audited_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                               retrieval_mode='existing_local' if offline else 'fresh_download',
                               protocol_version='0.1.0'))
    series, samples = parse_soft(paths['soft'])
    manifest = []
    for gsm, sample in samples.items():
        def one(key):
            values = sample[key]
            if len(values) != 1:
                raise ValueError('Missing/ambiguous ' + key + ' for ' + gsm)
            return values[0]
        title = one('!Sample_title')
        traits = dict(v.split(': ', 1) for v in sample['!Sample_characteristics_ch1'])
        treatment = traits.get('treatment')
        condition = {'vehicle': 'control', '10mg/kg oxaliplatin': 'oxaliplatin'}.get(treatment)
        if (condition is None or one('!Sample_organism_ch1') != 'Mus musculus'
                or traits.get('tissue') != 'Dorsal Root Ganglia (DRG)'):
            raise ValueError('Unexpected condition/species/tissue: ' + gsm)
        manifest.append(dict(sample_id=title, GEO_sample_id=gsm, condition=condition,
                             animal_or_pool_id='NA', biological_unit='independent_RNA_sample_reported',
                             species='Mus musculus', tissue=traits['tissue'], strain=traits.get('strain', 'NA'),
                             dose='10 mg/kg' if condition == 'oxaliplatin' else '0 mg/kg',
                             route='intraperitoneal', collection_time='after_8_weeks_treatment',
                             sex='male', age='NA', batch='NA',
                             inclusion_status='held', reason='R_sample_QC_and_metadata_review_pending',
                             study_level_source='GEO_series_summary_and_overall_design',
                             SRA_relation=';'.join(v for v in sample['!Sample_relation'] if v.startswith('SRA:')),
                             protocol_version='0.1.0'))
    titles = [s['sample_id'] for s in manifest]
    if len(titles) != len(set(titles)) or collections.Counter(s['condition'] for s in manifest) != {'control': 5, 'oxaliplatin': 5}:
        raise ValueError('Unexpected sample design')
    columns, genes, annotations, duplicates = read_counts(paths['counts'], titles)
    manifest.sort(key=lambda s: columns.index(s['sample_id']))
    table(out / 'sample_manifest.csv', manifest, list(manifest[0]))
    table(out / 'source_manifest.csv', provenance, list(provenance[0]))
    table(out / 'source_gene_annotations.csv', annotations, list(annotations[0]))
    table(out / 'duplicate_gene_audit.csv', duplicates,
          ['stable_gene_id', 'retained_source_row', 'duplicate_source_row', 'action'])
    table(out / 'count_matrix.tsv',
          (dict(stable_gene_id=k, **dict(zip(columns, g['values']))) for k, g in genes.items()),
          ['stable_gene_id'] + columns, '\t')
    table(out / 'gene_annotations.csv',
          (dict(stable_gene_id=k, symbol_candidates=';'.join(sorted(g['symbols'])),
                entrez_candidates=';'.join(sorted(g['entrez'])), biotype=g['biotype']) for k, g in genes.items()),
          ['stable_gene_id', 'symbol_candidates', 'entrez_candidates', 'biotype'])
    audit_fields = ['step_id', 'object_level', 'object_id', 'status', 'reason_code',
                    'reason_text', 'input_source', 'input_checksum', 'protocol_version']
    object_rows = [dict(step_id='01', object_level='sample', object_id=s['GEO_sample_id'],
                        status='held', reason_code='QC_PENDING', reason_text=s['reason'],
                        input_source=provenance[0]['source_url'], input_checksum=provenance[0]['sha256'],
                        protocol_version='0.1.0') for s in manifest]
    object_rows.extend(dict(step_id='01', object_level='gene', object_id=k, status='pass',
                            reason_code='COUNT_INTEGRITY_PASS',
                            reason_text='Unique stable ID with validated counts; not DE/sample-QC approval',
                            input_source=provenance[1]['source_url'], input_checksum=provenance[1]['sha256'],
                            protocol_version='0.1.0') for k in genes)
    table(out / 'step_audit.csv', object_rows, audit_fields)
    summary = dict(status='DATA_READY_QC_PENDING', accession='GSE286387', protocol_version='0.1.0',
                   input_rows=len(annotations), unique_genes=len(genes), n_samples=len(columns),
                   duplicate_gene_ids=len({r['stable_gene_id'] for r in duplicates}),
                   redundant_rows_removed=len(duplicates),
                   protein_coding_genes=sum(g['biotype'] == 'protein_coding' for g in genes.values()),
                   all_zero_genes=sum(max(g['values']) == 0 for g in genes.values()),
                   genes_max_below100=sum(max(g['values']) < 100 for g in genes.values()),
                   genes_passing_planned_filter=sum(sum(v >= 10 for v in g['values']) >= 5 for g in genes.values()),
                   library_totals={c: sum(g['values'][i] for g in genes.values()) for i, c in enumerate(columns)},
                   sources=provenance, python=sys.version, platform=platform.platform(),
                   remaining=['R sample-level QC before inclusion freeze',
                              'animal/pooling identities, age and batch not specified in GEO sample metadata',
                              'reference is Ensembl release 76; version-aware later mapping required',
                              'GEO mentions max-count filtering; export includes genes below 100, so stage of original filtering is unresolved'])
    (out / 'audit_summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
    report = '# Step 01 data audit\n\n' + '\n'.join(f'- {k}: {summary[k]}' for k in
                ['status', 'input_rows', 'unique_genes', 'n_samples', 'duplicate_gene_ids',
                 'redundant_rows_removed', 'protein_coding_genes', 'all_zero_genes', 'genes_passing_planned_filter'])
    report += '\n\nDuplicate stable IDs are collapsed only when all counts and biotype agree. Counts are never summed or rounded. Original annotation alternatives remain in source_gene_annotations.csv. No sample is excluded or approved before R QC. No DEGs are reported.\n\n'
    report += '\n'.join('- ' + item for item in summary['remaining']) + '\n'
    (out / 'qc_report.md').write_text(report, encoding='utf-8')
    print(report)
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-dir', type=Path, default=ROOT / 'data/raw/GSE286387')
    parser.add_argument('--out-dir', type=Path, default=ROOT / 'results/01_dataset_audit')
    parser.add_argument('--offline', action='store_true', help='Audit previously downloaded files; do not access network')
    args = parser.parse_args()
    try:
        run(args.raw_dir, args.out_dir, args.offline)
    except Exception as exc:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        for name in ['count_matrix.tsv', 'audit_summary.json', 'qc_report.md']:
            (args.out_dir / name).unlink(missing_ok=True)
        (args.out_dir / 'FAILURE.txt').write_text(str(exc) + '\n', encoding='utf-8')
        sys.exit('AUDIT FAILED: ' + str(exc))
    (args.out_dir / 'FAILURE.txt').unlink(missing_ok=True)
