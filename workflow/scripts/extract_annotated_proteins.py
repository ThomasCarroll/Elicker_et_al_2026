#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def parse_attrs(field):
    attrs = {}
    for part in field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k] = v
    return attrs


def load_gff(path: Path):
    mrnas = []
    cds_by_parent = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, _, feature, start, end, _, strand, _, attrs = parts
            attrs_map = parse_attrs(attrs)
            if feature == "mRNA" and "ID" in attrs_map:
                mrnas.append(
                    {
                        "id": attrs_map["ID"],
                        "seqid": seqid,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                    }
                )
            elif feature == "CDS" and "Parent" in attrs_map:
                cds_by_parent[attrs_map["Parent"]].append(line)
    return mrnas, cds_by_parent


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output-fasta", required=True)
    ap.add_argument("--output-bed", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--annotation-suffix", required=True)
    ap.add_argument("--genome-suffix", required=True)
    args = ap.parse_args()

    gff_cache = {}
    seen = set()

    Path(args.output_fasta).parent.mkdir(parents=True, exist_ok=True)
    with open(args.input) as in_handle, open(args.output_fasta, "w") as fasta_out, open(args.output_bed, "w") as bed_out:
        reader = csv.DictReader(in_handle, delimiter="\t")
        bed_out.write("contig\tstart\tend\tspecies\thit_id\treason\tstrand\n")

        for row in reader:
            species = row["species"]
            if species not in gff_cache:
                gff_cache[species] = load_gff(Path(args.genome_dir) / f"{species}{args.annotation_suffix}")
            mrnas, cds_by_parent = gff_cache[species]

            contig = row["chromosome"]
            start = int(row["sstart"])
            end = int(row["send"])
            low, high = sorted((start, end))

            overlapping = [
                mrna for mrna in mrnas
                if mrna["seqid"] == contig and mrna["start"] <= high and mrna["end"] >= low
            ]
            if not overlapping:
                bed_out.write(f"{contig}\t{low}\t{high}\t{species}\t{row['hit_id']}\tNo_mRNA_found\t{row['strand']}\n")
                continue

            best = None
            best_len = -1
            for mrna in overlapping:
                cds_len = 0
                for cds in cds_by_parent.get(mrna["id"], []):
                    parts = cds.split("\t")
                    cds_len += int(parts[4]) - int(parts[3]) + 1
                if cds_len > best_len:
                    best = mrna
                    best_len = cds_len

            if best is None or not cds_by_parent.get(best["id"]):
                bed_out.write(f"{contig}\t{low}\t{high}\t{species}\t{row['hit_id']}\tNo_CDS_found\t{row['strand']}\n")
                continue

            full_id = f"{species}_{best['id']}"
            if full_id in seen:
                continue
            seen.add(full_id)

            genome_fasta = Path(args.genome_dir) / f"{species}{args.genome_suffix}"
            with tempfile.TemporaryDirectory() as tmpdir:
                tmpdir = Path(tmpdir)
                gff_path = tmpdir / "cds.gff"
                cds_fasta = tmpdir / "cds.fa"
                protein_fasta = tmpdir / "protein.fa"
                with open(gff_path, "w") as gff_handle:
                    for cds_line in cds_by_parent[best["id"]]:
                        gff_handle.write(cds_line)
                subprocess.run(["gffread", "-g", str(genome_fasta), "-x", str(cds_fasta), str(gff_path)], check=True)
                subprocess.run(["transeq", "-sequence", str(cds_fasta), "-outseq", str(protein_fasta), "-frame", "1", "-clean"], check=True)

                records = list(fasta_records(protein_fasta))
                if not records:
                    bed_out.write(f"{contig}\t{low}\t{high}\t{species}\t{row['hit_id']}\tTranslation_failed\t{row['strand']}\n")
                    continue
                fasta_out.write(f">{full_id}\n{records[0][1]}\n")


if __name__ == "__main__":
    main()
