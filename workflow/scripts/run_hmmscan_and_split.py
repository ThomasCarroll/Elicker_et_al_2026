#!/usr/bin/env python3

import argparse
import csv
import shutil
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


def fasta_has_records(path: Path) -> bool:
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                return True
    return False


def read_species_list(path: Path):
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def write_empty_outputs(args):
    Path(args.domtbl).write_text("")
    with open(args.classification, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["query", "species", "gene", "domains_detected", "domain_counts", "classification"])
    with open(args.species_summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["species", "FULL", "PARTIAL", "IGNORE"])
        for species in read_species_list(Path(args.species_list)):
            writer.writerow([species, 0, 0, 0])
    Path(args.full).write_text("")
    Path(args.partial).write_text("")
    Path(args.ignore).write_text("")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--input", required=True)
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--domtbl", required=True)
    ap.add_argument("--classification", required=True)
    ap.add_argument("--species-summary", required=True)
    ap.add_argument("--species-list", required=True)
    ap.add_argument("--classifier", required=True)
    ap.add_argument("--full", required=True)
    ap.add_argument("--partial", required=True)
    ap.add_argument("--ignore", required=True)
    args = ap.parse_args()

    Path(args.domtbl).parent.mkdir(parents=True, exist_ok=True)
    Path(args.classification).parent.mkdir(parents=True, exist_ok=True)
    Path(args.species_summary).parent.mkdir(parents=True, exist_ok=True)

    if not fasta_has_records(Path(args.input)):
        write_empty_outputs(args)
        return

    with open(Path(args.domtbl).parent / f"{args.gene}_hmmscan.log", "w") as log_handle:
        subprocess.run(
            ["hmmscan", "--cpu", "8", "--domtblout", args.domtbl, args.pfam, args.input],
            check=True,
            stdout=log_handle,
        )

    subprocess.run(
        ["python3", args.classifier, args.gene, args.domtbl, args.classification, args.species_list],
        check=True,
    )

    generated_species_summary = Path(f"{args.gene}_species_summary.tsv")
    Path(args.species_summary).parent.mkdir(parents=True, exist_ok=True)
    if generated_species_summary.exists():
        shutil.move(str(generated_species_summary), args.species_summary)

    class_map = {}
    with open(args.classification) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            class_map[row["query"]] = row["classification"].strip()

    with open(args.full, "w") as full, open(args.partial, "w") as partial, open(args.ignore, "w") as ignore:
        for name, seq in fasta_records(Path(args.input)):
            klass = class_map.get(name, "IGNORE")
            out = {"FULL": full, "PARTIAL": partial}.get(klass, ignore)
            out.write(f">{name}\n{seq}\n")


if __name__ == "__main__":
    main()
