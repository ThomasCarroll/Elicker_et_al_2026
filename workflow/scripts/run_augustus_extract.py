#!/usr/bin/env python3

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def clean_augustus_output(path: Path, out_fasta: Path):
    seq = ""
    header = ""
    in_seq = False
    with open(path) as handle, open(out_fasta, "w") as out:
        for line in handle:
            if line.startswith("# ----- prediction on sequence number"):
                if seq and header:
                    out.write(f">{header}\n{seq}\n")
                seq = ""
                header = ""
                in_seq = False
            elif "name =" in line:
                marker = "name = "
                header = line.split(marker, 1)[1].strip()
            elif line.startswith("# protein sequence = ["):
                in_seq = True
                seq = line.replace("# protein sequence = [", "").split("]")[0].strip()
            elif in_seq and line.startswith("#"):
                clean = line.lstrip("# ").rstrip("\n")
                if "]" in clean:
                    seq += clean.split("]")[0]
                    in_seq = False
                else:
                    seq += clean
        if seq and header:
            out.write(f">{header}\n{seq}\n")


def ensure_fai(genome: Path):
    fai = Path(f"{genome}.fai")
    if not fai.exists() or fai.stat().st_mtime < genome.stat().st_mtime:
        subprocess.run(["samtools", "faidx", str(genome)], check=True)
    return fai


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-bed", required=True)
    ap.add_argument("--output-fasta", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", required=True)
    ap.add_argument("--augustus-species", required=True)
    ap.add_argument("--extend-bp", type=int, required=True)
    args = ap.parse_args()

    outdir = Path(args.output_fasta).parent
    outdir.mkdir(parents=True, exist_ok=True)
    sequences = outdir / "augustus_regions.fa"

    with open(args.input_bed) as in_handle, open(sequences, "w") as seq_out:
        reader = csv.DictReader(in_handle, delimiter="\t")
        for row in reader:
            species = row["species"]
            chrom = row["contig"]
            start = int(row["start"])
            end = int(row["end"])
            strand = row["strand"]
            genome = Path(args.genome_dir) / f"{species}{args.genome_suffix}"
            fai_path = ensure_fai(genome)
            contig_len = None
            with open(fai_path) as fai:
                for line in fai:
                    parts = line.split("\t")
                    if parts[0] == chrom:
                        contig_len = int(parts[1])
                        break
            if contig_len is None:
                continue
            ext_start = max(0, min(start, end) - args.extend_bp)
            ext_end = min(contig_len, max(start, end) + args.extend_bp)
            with tempfile.NamedTemporaryFile("w", delete=False) as bed_handle:
                bed_handle.write(f"{chrom}\t{ext_start}\t{ext_end}\t{species}|{chrom}:{ext_start}-{ext_end}({strand})\t0\t{strand}\n")
                bed_path = bed_handle.name
            subprocess.run(["bedtools", "getfasta", "-fi", str(genome), "-bed", bed_path, "-s", "-name"], check=True, stdout=seq_out)
            Path(bed_path).unlink(missing_ok=True)

    augustus_out = outdir / "augustus_output.txt"
    with open(augustus_out, "w") as handle:
        subprocess.run(["augustus", f"--species={args.augustus_species}", str(sequences), "--gff3=on", "--protein=on"], check=True, stdout=handle)
    clean_augustus_output(augustus_out, Path(args.output_fasta))


if __name__ == "__main__":
    main()
