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


def parse_attrs(field):
    attrs = {}
    for part in field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def read_threshold(path: Path) -> int:
    with open(path) as handle:
        for line in handle:
            if line.startswith("Overall longest intron threshold"):
                return int(line.rstrip("\n").split("\t")[1])
    return 0


def ensure_fai(genome: Path):
    fai = Path(f"{genome}.fai")
    if not fai.exists() or fai.stat().st_mtime < genome.stat().st_mtime:
        subprocess.run(["samtools", "faidx", str(genome)], check=True)
    return fai


def contig_lengths(path: Path):
    values = {}
    with open(path) as handle:
        for line in handle:
            contig, length, *_ = line.rstrip("\n").split("\t")
            values[contig] = int(length)
    return values


def gene_coords(gff: Path, transcript_id: str):
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                continue
            attrs = parse_attrs(parts[8])
            if attrs.get("ID") == transcript_id:
                return parts[0], int(parts[3]), int(parts[4])
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["annotated", "predicted"], required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--threshold-file", required=True)
    ap.add_argument("--likely-full", required=True)
    ap.add_argument("--remaining", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--annotation-suffix", required=True)
    ap.add_argument("--genome-suffix", required=True)
    args = ap.parse_args()

    threshold = read_threshold(Path(args.threshold_file))
    Path(args.summary).parent.mkdir(parents=True, exist_ok=True)

    with open(args.likely_full, "w") as likely, open(args.remaining, "w") as remaining, open(args.summary, "w") as summary:
        summary.write("species\tgene_id\tcontig\tstart\tend\tcontig_len\tnear_end\n")
        for header, seq in fasta_records(Path(args.input)):
            if args.mode == "annotated":
                parts = header.split("_")
                gene_id = "_".join(parts[-2:])
                species = "_".join(parts[:-2])
                coords = gene_coords(Path(args.genome_dir) / f"{species}{args.annotation_suffix}", gene_id)
                if coords is None:
                    continue
                contig, start, end = coords
            else:
                species, remainder = header.split("|", 1)
                contig, coords = remainder.split(":")
                region, _strand = coords.split("(")
                start, end = [int(x) for x in region.rstrip(")").split("-")]
                gene_id = header

            genome = Path(args.genome_dir) / f"{species}{args.genome_suffix}"
            fai = ensure_fai(genome)
            lengths = contig_lengths(fai)
            contig_len = lengths.get(contig)
            if contig_len is None:
                continue

            near_start = start <= threshold
            near_end = (contig_len - end) <= threshold
            out = likely if (near_start or near_end) else remaining
            out.write(f">{header}\n{seq}\n")
            summary.write(f"{species}\t{gene_id}\t{contig}\t{start}\t{end}\t{contig_len}\t{'YES' if (near_start or near_end) else 'NO'}\n")


if __name__ == "__main__":
    main()
