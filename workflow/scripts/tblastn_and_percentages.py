#!/usr/bin/env python3

import argparse
import csv
import subprocess
from collections import defaultdict
from pathlib import Path


HEADER = [
    "query", "species", "chromosome", "qstart", "qend", "sstart", "send",
    "length", "evalue", "bitscore", "pident", "qcovs", "sseq", "strand"
]


def read_fasta_length(path: Path) -> int:
    total = 0
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def ensure_blast_db(genome: Path, db_prefix: Path) -> None:
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["makeblastdb", "-in", str(genome), "-dbtype", "nucl", "-out", str(db_prefix)],
        check=True,
        capture_output=True,
        text=True,
    )


def run_tblastn(query: Path, db_prefix: Path):
    cmd = [
        "tblastn",
        "-query", str(query),
        "-db", str(db_prefix),
        "-outfmt",
        "6 qseqid sseqid qstart qend sstart send length evalue bitscore pident qcovs sseq",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"tblastn failed for database {db_prefix}\n"
            f"command: {' '.join(cmd)}\n"
            f"stderr:\n{proc.stderr}"
        )
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 12:
            continue
        strand = "-" if int(fields[5]) < int(fields[4]) else "+"
        yield [
            fields[0], None, fields[1], int(fields[2]), int(fields[3]), int(fields[4]),
            int(fields[5]), int(fields[6]), fields[7], fields[8], float(fields[9]),
            float(fields[10]), fields[11], strand,
        ]


def write_tsv(path: Path, rows, header):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", required=True)
    ap.add_argument("--gene", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-glob", required=True)
    ap.add_argument("--combined", required=True)
    ap.add_argument("--sorted", required=True)
    ap.add_argument("--grouped", required=True)
    ap.add_argument("--identity", required=True)
    ap.add_argument("--coverage", required=True)
    args = ap.parse_args()

    query = Path(args.query)
    genome_dir = Path(args.genome_dir)
    combined_rows = []
    blast_db_dir = Path(args.combined).parent / "_blastdb"

    for genome in sorted(genome_dir.glob(args.genome_glob)):
        species = genome.name.replace("_genome.fasta", "")
        db_prefix = blast_db_dir / f"{species}_genome"
        ensure_blast_db(genome, db_prefix)
        for row in run_tblastn(query, db_prefix):
            row[1] = species
            combined_rows.append(row)

    write_tsv(Path(args.combined), combined_rows, HEADER)

    sorted_rows = sorted(combined_rows, key=lambda r: (r[1], r[2], r[5]))
    write_tsv(Path(args.sorted), sorted_rows, HEADER)

    grouped_header = HEADER[:2] + ["hit_id"] + HEADER[2:]
    grouped_rows = []
    hit_id = 0
    prev = None
    for row in sorted_rows:
        key = (row[1], row[2], row[13])
        qstart = row[3]
        if prev is None or key != prev[:3]:
            hit_id += 1
        else:
            prev_qstart = prev[3]
            if row[13] == "+" and qstart <= prev_qstart:
                hit_id += 1
            if row[13] == "-" and qstart >= prev_qstart:
                hit_id += 1
        grouped_rows.append(row[:2] + [hit_id] + row[2:])
        prev = (row[1], row[2], row[13], qstart)
    write_tsv(Path(args.grouped), grouped_rows, grouped_header)

    identity_header = grouped_header + ["hit_percent_identity"]
    by_hit = defaultdict(list)
    for row in grouped_rows:
        by_hit[row[2]].append(row)
    identity_rows = []
    for hit_rows in by_hit.values():
        total_len = sum(r[8] for r in hit_rows)
        weighted = sum(r[8] * r[11] for r in hit_rows) / total_len if total_len else 0
        for row in hit_rows:
            identity_rows.append(row + [round(weighted, 1)])
    identity_rows.sort(key=lambda r: (r[1], r[3], r[6]))
    write_tsv(Path(args.identity), identity_rows, identity_header)

    coverage_header = identity_header + ["hit_coverage_percent"]
    query_length = read_fasta_length(query)
    coverage_rows = []
    for hit, hit_rows in defaultdict(list, {k: [r for r in identity_rows if r[2] == k] for k in {r[2] for r in identity_rows}}).items():
        intervals = []
        for row in hit_rows:
            start, end = sorted((row[4], row[5]))
            intervals.append((start, end))
        intervals.sort()
        merged = []
        for start, end in intervals:
            if not merged or start > merged[-1][1]:
                merged.append([start, end])
            else:
                merged[-1][1] = max(merged[-1][1], end)
        covered = sum(end - start + 1 for start, end in merged)
        coverage = round((covered / query_length) * 100, 1) if query_length else 0.0
        for row in hit_rows:
            coverage_rows.append(row + [coverage])
    coverage_rows.sort(key=lambda r: (r[1], r[3], r[6]))
    write_tsv(Path(args.coverage), coverage_rows, coverage_header)


if __name__ == "__main__":
    main()
