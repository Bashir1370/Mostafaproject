#!/usr/bin/env python3
"""Inventory Cys in reviewed mouse canonical sequences; offline, no scoring."""
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import sys
import tempfile

VERSION = '1.0.0'
ALPHABET = set('ACDEFGHIKLMNPQRSTVWYOU')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fields):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def safe_path(project, relative):
    path = (project / relative).resolve()
    require(project in path.parents, 'Input path outside project: ' + relative)
    return path


def read_fasta(path):
    sequences, metadata = {}, {}
    accession = None
    for line in path.read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            match = re.fullmatch(r'>(\S+) taxid=10090 genes=(ENSMUSG\d{11}(?:,ENSMUSG\d{11})*)', line)
            require(match is not None, 'Invalid mouse FASTA header: ' + line)
            accession, genes = match.groups()
            require(accession not in sequences, 'Duplicate FASTA accession: ' + accession)
            gene_list = genes.split(',')
            require(len(gene_list) == len(set(gene_list)), 'Duplicate FASTA gene association')
            sequences[accession] = ''
            metadata[accession] = set(gene_list)
        else:
            require(accession is not None and bool(line) and set(line) <= ALPHABET,
                    'Empty, ambiguous or invalid FASTA sequence line')
            sequences[accession] += line
    require(bool(sequences) and all(sequences.values()), 'Missing/empty FASTA sequence')
    return sequences, metadata


def load_input(project):
    binding_path = project / 'config/GSE286387_step04_input.json'
    binding = json.loads(binding_path.read_text(encoding='utf-8'))
    require(binding['status'] == 'STEP03B_SNAPSHOT_ACCEPTED_FOR_STEP04', 'Mapping review not accepted')
    require(binding['protocol_version'] == '0.1.0' and binding['species_taxid'] == 10090,
            'Unexpected protocol/species')
    inputs = [binding_path]
    for relative, expected in binding['input_sha256'].items():
        path = safe_path(project, relative)
        require(path.is_file() and sha(path.read_bytes()) == expected, 'Frozen input checksum changed: ' + relative)
        inputs.append(path)
    for relative in (binding['mapping_file'], binding['fasta_file'], binding['discovery_file']):
        require(relative in binding['input_sha256'], 'Unbound primary input: ' + relative)
    config = json.loads((project / 'config/analysis_parameters.json').read_text())
    require(config['protocol_version'] == binding['protocol_version'] and
            config['sequence']['position_numbering'] == '1_based' and
            config['scope']['taxonomy_id'] == 10090, 'Protocol numbering/species changed')
    rows = read_csv(safe_path(project, binding['mapping_file']))
    discovery = read_csv(safe_path(project, binding['discovery_file']))
    genes = [r['stable_gene_id'] for r in rows]
    require(len(rows) == binding['input_genes'] and len(genes) == len(set(genes)) and
            set(genes) == {r['stable_gene_id'] for r in discovery}, 'Discovery/mapping universe mismatch')
    sequences, metadata = read_fasta(safe_path(project, binding['fasta_file']))
    associated = defaultdict(set)
    for row in rows:
        status = row['mapping_status']
        require(status in ('pass', 'held', 'unassessable'), 'Invalid mapping status')
        if status == 'pass':
            accession = row['protein_accession']
            require(accession in sequences and row['species_taxid'] == '10090' and
                    row['canonical_status'] == 'UniProtKB_displayed_sequence' and
                    row['reference_release'] == binding['uniprot_sequence_release'] and
                    row['protocol_version'] == binding['protocol_version'], 'Accepted identity mismatch')
            sequence = sequences[accession]
            require(len(sequence) == int(row['sequence_length']) and
                    sha(sequence.encode()) == row['sequence_checksum'], 'Sequence length/checksum mismatch: ' + accession)
            associated[accession].add(row['stable_gene_id'])
        else:
            require(not row['protein_accession'] and not row['sequence_checksum'], 'Unaccepted mapping claims sequence')
    require(dict(associated) == metadata and set(sequences) == set(associated), 'FASTA gene/accession associations mismatch')
    require(sum(r['mapping_status'] == 'pass' for r in rows) == binding['passed_genes'] and
            len(sequences) == binding['unique_accepted_proteins'], 'Accepted coverage mismatch')
    return binding, rows, sequences, associated, inputs


def inventory(sequences, associated, release, protocol):
    proteins, sites = [], []
    for accession, sequence in sorted(sequences.items()):
        positions = [i for i, residue in enumerate(sequence, 1) if residue == 'C']
        count = len(positions)
        common = dict(protein_accession=accession, stable_gene_ids=';'.join(sorted(associated[accession])),
                      sequence_length=len(sequence), sequence_checksum=sha(sequence.encode()),
                      n_cys_total=count, cys_percent=100.0 * count / len(sequence),
                      reference_release=release, protocol_version=protocol)
        proteins.append(dict(common, cys_positions=';'.join(map(str, positions)),
                             inventory_status='pass' if count else 'excluded',
                             reason='AT_LEAST_ONE_CYS' if count else 'NO_CYSTEINE_IN_ACCEPTED_SEQUENCE'))
        for position in positions:
            sites.append(dict(common, site_id=f'{accession}:C{position}', canonical_cys_position=position))
    return proteins, sites


def run(project):
    out = project / 'results/04_cysteine_inventory'
    out.mkdir(parents=True, exist_ok=True)
    (out / 'SUCCESS.txt').unlink(missing_ok=True)
    stage = Path(tempfile.mkdtemp(prefix='inventory_stage_', dir=out))
    try:
        print('[1/3] Verifying reviewed snapshot checksums and canonical FASTA...', flush=True)
        binding, rows, sequences, associated, inputs = load_input(project)
        print('[2/3] Counting Cys and recording 1-based canonical positions...', flush=True)
        proteins, sites = inventory(sequences, associated, binding['uniprot_sequence_release'], binding['protocol_version'])
        protein_fields = ['protein_accession', 'stable_gene_ids', 'sequence_length', 'sequence_checksum',
                          'n_cys_total', 'cys_percent', 'reference_release', 'protocol_version',
                          'cys_positions', 'inventory_status', 'reason']
        site_fields = ['site_id', 'protein_accession', 'stable_gene_ids', 'canonical_cys_position',
                       'n_cys_total', 'cys_percent', 'sequence_length', 'sequence_checksum',
                       'reference_release', 'protocol_version']
        write_csv(stage / 'protein_inventory.csv', proteins, protein_fields)
        write_csv(stage / 'cysteine_inventory.csv', sites, site_fields)
        indexed = {r['protein_accession']: r for r in proteins}
        mapping_sha = binding['input_sha256'][binding['mapping_file']]
        fasta_sha = binding['input_sha256'][binding['fasta_file']]
        audit_fields = ['step_id', 'object_level', 'object_id', 'status', 'reason_code', 'reason_text',
                        'input_source', 'input_checksum', 'protocol_version']
        audit = [dict(step_id='04', object_level='protein', object_id=r['protein_accession'],
                      status=r['inventory_status'], reason_code=r['reason'], reason_text=r['reason'],
                      input_source=binding['fasta_file'], input_checksum=fasta_sha,
                      protocol_version=binding['protocol_version']) for r in proteins]
        write_csv(stage / 'step_audit.csv', audit, audit_fields)
        gene_audit = []
        for row in rows:
            if row['mapping_status'] == 'pass':
                protein = indexed[row['protein_accession']]
                status, reason = protein['inventory_status'], protein['reason']
            else:
                status, reason = row['mapping_status'], 'UPSTREAM_MAPPING_' + row['reason']
            gene_audit.append(dict(step_id='04', object_level='gene', object_id=row['stable_gene_id'],
                                   status=status, reason_code=reason, reason_text=reason,
                                   input_source=binding['mapping_file'], input_checksum=mapping_sha,
                                   protocol_version=binding['protocol_version']))
        write_csv(stage / 'gene_step_audit.csv', gene_audit, audit_fields)
        with (stage / 'cysteine_positive_sequences.fasta').open('w', encoding='utf-8') as handle:
            for accession, sequence in sorted(sequences.items()):
                if indexed[accession]['inventory_status'] != 'pass':
                    continue
                handle.write('>' + accession + ' taxid=10090 genes=' + ','.join(sorted(associated[accession])) + '\n')
                handle.write('\n'.join(sequence[i:i + 60] for i in range(0, len(sequence), 60)) + '\n')
        inputs.append(Path(__file__).resolve())
        write_csv(stage / 'input_checksums.csv',
                  [dict(input=str(p.relative_to(project)), sha256=sha(p.read_bytes())) for p in inputs],
                  ['input', 'sha256'])
        summary = dict(status='STEP04_GENERATED_REVIEW_PENDING', protocol_version=binding['protocol_version'],
                       implementation_version=VERSION, input_discovery_genes=len(rows),
                       input_mapped_genes=binding['passed_genes'], input_unique_proteins=len(proteins),
                       protein_status_counts=dict(Counter(r['inventory_status'] for r in proteins)),
                       gene_status_counts=dict(Counter(r['status'] for r in gene_audit)),
                       total_unique_cysteine_sites=len(sites),
                       proteins_with_one_cys=sum(r['n_cys_total'] == 1 for r in proteins),
                       proteins_with_multiple_cys=sum(r['n_cys_total'] > 1 for r in proteins),
                       uniprot_sequence_release=binding['uniprot_sequence_release'],
                       mapping_sha256=mapping_sha, fasta_sha256=fasta_sha,
                       python=platform.python_version(), platform=platform.platform(),
                       generated_utc=datetime.now(timezone.utc).isoformat())
        (stage / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n', encoding='utf-8')
        report = ['# Step 04 — canonical cysteine inventory', '', 'Status: ' + summary['status'], '',
                  f'- Discovery genes: {len(rows)}; mapped genes: {binding["passed_genes"]}',
                  f'- Unique input proteins: {len(proteins)}',
                  f'- Proteins with Cys: {summary["protein_status_counts"].get("pass", 0)}',
                  f'- Proteins without Cys: {summary["protein_status_counts"].get("excluded", 0)}',
                  f'- Proteins with exactly one Cys: {summary["proteins_with_one_cys"]}',
                  f'- Proteins with multiple Cys: {summary["proteins_with_multiple_cys"]}',
                  f'- Unique Cys sites: {len(sites)}', '',
                  'Positions are 1-based in the unchanged UniProtKB displayed sequence, including precursor regions.',
                  'U (selenocysteine) is retained in sequence but is not counted as Cys (C).',
                  'Each accession/site is counted once, even when associated with multiple stable genes.',
                  'Cys count and percentage are descriptive and do not contribute to the structural score.',
                  'No-Cys proteins are excluded only from this Cys-based analysis, not declared oxidation-resistant.',
                  'Upstream held/unassessable genes remain in gene_step_audit.csv and were not assessed for Cys.',
                  'Chemical state, mature-chain processing, structural coverage, SASA and pKa have not been assessed.',
                  'Review workstation outputs before binding Step 05 input.']
        (stage / 'inventory_report.md').write_text('\n'.join(report) + '\n', encoding='utf-8')
        for path in stage.iterdir():
            shutil.copyfile(path, out / path.name)
        (out / 'FAILURE.txt').unlink(missing_ok=True)
        (out / 'SUCCESS.txt').write_text(summary['status'] + '\n', encoding='utf-8')
        print('[3/3] Outputs: ' + str(out), flush=True)
        print('\n'.join(report[:11]))
        return summary
    except Exception as error:
        (out / 'FAILURE.txt').write_text(str(error) + '\n', encoding='utf-8')
        raise
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
        print('STEP04 FAILED: ' + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
