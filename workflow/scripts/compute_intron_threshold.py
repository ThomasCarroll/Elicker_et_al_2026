#!/usr/bin/env python3

import argparse
from pathlib import Path


def quantile(values, q):
    if not values:
        return 0
    pos = (len(values) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(values) - 1)
    frac = pos - lo
    return round(values[lo] + (values[hi] - values[lo]) * frac)


def fasta_headers(path: Path):
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                yield line[1:].strip()


def parse_attrs(field):
    attrs = {}
    for part in field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def longest_intron(gff: Path, transcript: str) -> int:
    longest = 0
    with open(gff) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "intron":
                continue
            attrs = parse_attrs(parts[8])
            if attrs.get("Parent") == transcript:
                longest = max(longest, int(parts[4]) - int(parts[3]) + 1)
    return longest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--annotation-suffix", required=True)
    args = ap.parse_args()

    lengths = []
    rows = []
    for header in fasta_headers(Path(args.input)):
        parts = header.split("_")
        transcript = "_".join(parts[-2:])
        species = "_".join(parts[:-2])
        gff = Path(args.genome_dir) / f"{species}{args.annotation_suffix}"
        value = longest_intron(gff, transcript) if gff.exists() else 0
        lengths.append(value)
        rows.append((transcript, value))

    lengths.sort()
    q1 = quantile(lengths, 0.25)
    q3 = quantile(lengths, 0.75)
    iqr = q3 - q1
    threshold = q3 + 3 * iqr

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as handle:
        for transcript, value in rows:
            handle.write(f"{transcript}\t{value}\n")
        handle.write(f"\nOverall longest intron threshold\t{threshold}\n")


if __name__ == "__main__":
    main()
