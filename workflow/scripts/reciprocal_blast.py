#!/usr/bin/env python3

import argparse
import subprocess
from pathlib import Path


def fasta_records(path: Path):
    name, seq = None, []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def parse_gene_name(attrs: str):
    for item in attrs.split(";"):
        if item.startswith("gene=") or item.startswith("Gene="):
            return item.split("=", 1)[1]
    return ""


def overlapping_gene(gff: Path, contig: str, start: int, end: int):
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[0] != contig or parts[2] != "gene":
                continue
            if int(parts[3]) <= end and int(parts[4]) >= start:
                return parse_gene_name(parts[8])
    return ""


def ensure_blast_db(reference_fasta: Path, db_prefix: Path):
    db_prefix.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["makeblastdb", "-in", str(reference_fasta), "-dbtype", "nucl", "-out", str(db_prefix)],
        check=True,
        capture_output=True,
        text=True,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--confirmed", required=True)
    ap.add_argument("--eliminated", required=True)
    ap.add_argument("--reference-fasta", required=True)
    ap.add_argument("--blast-db", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--gene-match", required=True)
    args = ap.parse_args()

    Path(args.confirmed).parent.mkdir(parents=True, exist_ok=True)
    db_prefix = Path(args.confirmed).parent / "_blastdb" / Path(args.blast_db).name
    ensure_blast_db(Path(args.reference_fasta), db_prefix)
    with open(args.confirmed, "w") as confirmed, open(args.eliminated, "w") as eliminated:
        for name, seq in fasta_records(Path(args.input)):
            tmp = Path(args.confirmed).parent / f".{name}.query.fa"
            tmp.write_text(f">{name}\n{seq}\n")
            cmd = [
                "tblastn", "-query", str(tmp), "-db", str(db_prefix), "-evalue", "1e-5",
                "-outfmt", "6 sseqid sstart send", "-max_target_seqs", "1",
            ]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            tmp.unlink(missing_ok=True)
            if proc.returncode != 0:
                raise RuntimeError(
                    f"tblastn failed for reciprocal search\n"
                    f"command: {' '.join(cmd)}\n"
                    f"stderr:\n{proc.stderr}"
                )
            if not proc.stdout.strip():
                eliminated.write(f">{name}\n{seq}\n")
                continue
            contig, start, end = proc.stdout.splitlines()[0].split("\t")
            low, high = sorted((int(start), int(end)))
            gene_name = overlapping_gene(Path(args.gff), contig, low, high)
            if args.gene_match.lower() in gene_name.lower():
                confirmed.write(f">{name}\n{seq}\n")
            else:
                eliminated.write(f">{name}\n{seq}\n")


if __name__ == "__main__":
    main()
