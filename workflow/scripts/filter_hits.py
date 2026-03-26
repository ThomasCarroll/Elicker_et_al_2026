#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--identity-threshold", type=float, required=True)
    ap.add_argument("--coverage-threshold", type=float, required=True)
    args = ap.parse_args()

    rows = []
    with open(args.input) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames
        for row in reader:
            if float(row["hit_percent_identity"]) > args.identity_threshold and float(row["hit_coverage_percent"]) > args.coverage_threshold:
                rows.append(row)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
