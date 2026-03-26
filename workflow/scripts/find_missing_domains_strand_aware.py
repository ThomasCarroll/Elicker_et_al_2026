#!/usr/bin/env python3

import argparse
import csv
import subprocess
import tempfile
from collections import Counter
from pathlib import Path


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(table)[::-1]


def translate_three_frames(seq: str, strand: str):
    seq = seq.upper()
    if strand == "-":
        seq = revcomp(seq)
    records = []
    for frame in range(3):
        peptide = []
        for i in range(frame, len(seq) - 2, 3):
            peptide.append(CODON_TABLE.get(seq[i:i + 3], "X"))
        records.append((frame + 1, "".join(peptide).rstrip("*")))
    return records


def parse_counts(field: str):
    counts = Counter()
    if field and field != ".":
        for part in field.split(","):
            part = part.strip()
            if not part:
                continue
            if "(" in part:
                name, count = part.split("(", 1)
                counts[name] += int(count.rstrip(")"))
            else:
                counts[part] += 1
    return counts


def rules_for_gene(gene: str):
    required = {}
    either_or = []
    mixed = []
    if gene in ["AGO1", "AGO2"]:
        required = {"Piwi": 1, "ArgoMid": 1, "ArgoN": 1, "PAZ": 1, "ArgoL1": 1, "ArgoL2": 1}
    elif gene in ["AGO3", "PIWI"]:
        required = {"Piwi": 1, "PAZ": 1, "Piwi_N": 1}
    elif gene == "AUB":
        required = {"Piwi": 1, "ArgoN": 1, "PAZ": 1}
    elif gene == "PASHA":
        required = {"dsrm": 1}
    elif gene in ["LOQS", "R2D2"]:
        required = {"dsrm": 2}
    elif gene == "DROSHA":
        required = {"dsrm": 1}
        mixed = [{"domains": ["Ribonuclease_3", "Ribonuclease_3_3", "Ribonucleas_3_3"], "count": 2}]
    elif gene == "DICER1":
        required = {"Dicer_PBD": 1, "Helicase_C": 1, "Dicer_dimer": 1, "Dicer_platform": 1, "PAZ": 1, "Dicer_dsRBD": 1}
        mixed = [{"domains": ["Ribonuclease_3", "Ribonuclease_3_3", "Ribonucleas_3_3"], "count": 2}]
    elif gene == "DICER2":
        required = {"Dicer_dimer": 1, "Dicer_platform": 1, "PAZ": 1, "Helicase_C": 1, "Dicer_dsRBD": 1}
        either_or = [["DEAD", "ResIII"]]
        mixed = [{"domains": ["Ribonuclease_3", "Ribonuclease_3_3", "Ribonucleas_3_3"], "count": 2}]
    return required, either_or, mixed


def load_contig_length(fai: Path, contig: str):
    with open(fai) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == contig:
                return int(parts[1])
    return None


def ensure_fai(genome: Path):
    fai = Path(f"{genome}.fai")
    if not fai.exists() or fai.stat().st_mtime < genome.stat().st_mtime:
        subprocess.run(["samtools", "faidx", str(genome)], check=True)
    return fai


def extract_region(genome: Path, contig: str, start: int, end: int) -> str:
    proc = subprocess.run(["samtools", "faidx", str(genome), f"{contig}:{start}-{end}"], check=True, capture_output=True, text=True)
    return "".join(line.strip() for line in proc.stdout.splitlines() if not line.startswith(">"))


def hmmer_counts(pfam: Path, fasta_path: Path):
    with tempfile.NamedTemporaryFile("w+", delete=False) as domtbl:
        subprocess.run(["hmmscan", "--domtblout", domtbl.name, str(pfam), str(fasta_path)], check=True, stdout=subprocess.DEVNULL)
        counts = Counter()
        with open(domtbl.name) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) >= 23 and float(fields[6]) < 1e-5 and float(fields[11]) < 1e-5:
                    counts[fields[0]] += 1
    Path(domtbl.name).unlink(missing_ok=True)
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--input-table", required=True)
    ap.add_argument("--output-table", required=True)
    ap.add_argument("--pfam", required=True)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--genome-suffix", required=True)
    ap.add_argument("--strand-column", default="STRAND")
    args = ap.parse_args()

    required, either_or, mixed = rules_for_gene(args.gene)
    arch_domains = set(required)
    for group in either_or:
        arch_domains.update(group)
    for group in mixed:
        arch_domains.update(group["domains"])

    Path(args.output_table).parent.mkdir(parents=True, exist_ok=True)
    with open(args.input_table) as in_handle, open(args.output_table, "w", newline="") as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = reader.fieldnames + ["FOUND_DOMAINS", "RECLASSIFICATION"]
        writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            original = row.get("CLASSIFICATION", ".")
            if original != "PARTIAL":
                row["FOUND_DOMAINS"] = "."
                row["RECLASSIFICATION"] = "FULL" if original == "FULL" else original
                writer.writerow(row)
                continue

            species = row["SPECIES"]
            contig = row["CHROMOSOME"]
            start = int(row["START"])
            end = int(row["END"])
            strand = row.get(args.strand_column, row.get("STRAND", "+"))
            genome = Path(args.genome_dir) / f"{species}{args.genome_suffix}"
            contig_len = load_contig_length(ensure_fai(genome), contig)
            if contig_len is None:
                row["FOUND_DOMAINS"] = "."
                row["RECLASSIFICATION"] = "UNKNOWN_CONTIG"
                writer.writerow(row)
                continue

            ext_start = max(1, min(start, end) - 5000)
            ext_end = min(contig_len, max(start, end) + 5000)
            seq = extract_region(genome, contig, ext_start, ext_end)
            with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fa") as tmp_fa:
                for frame, peptide in translate_three_frames(seq, strand):
                    tmp_fa.write(f">{species}_{row['HIT_ID']}_frame{frame}\n{peptide}\n")
                tmp_path = Path(tmp_fa.name)

            counts = hmmer_counts(Path(args.pfam), tmp_path)
            tmp_path.unlink(missing_ok=True)
            relevant = Counter({k: v for k, v in counts.items() if k in arch_domains})
            annotated = parse_counts(row.get("DOMAINS", ""))

            contributed = False
            satisfied = True
            for dom, req in required.items():
                if annotated[dom] + relevant[dom] < req:
                    satisfied = False
                if min(req, annotated[dom] + relevant[dom]) > min(req, annotated[dom]):
                    contributed = True
            for group in mixed:
                need = group["count"]
                have_annotated = sum(annotated[d] for d in group["domains"])
                have_new = sum(relevant[d] for d in group["domains"])
                if have_annotated + have_new < need:
                    satisfied = False
                if min(need, have_annotated + have_new) > min(need, have_annotated):
                    contributed = True
            for group in either_or:
                old = any(annotated[d] > 0 for d in group)
                new = any(relevant[d] > 0 for d in group)
                if not (old or new):
                    satisfied = False
                if not old and new:
                    contributed = True

            row["FOUND_DOMAINS"] = ",".join(f"{k}({relevant[k]})" for k in sorted(relevant)) or "."
            if satisfied:
                row["RECLASSIFICATION"] = "LIKELY_FULL_DOMAINS_FOUND"
            elif contributed:
                row["RECLASSIFICATION"] = "PARTIAL_DOMAINS_FOUND"
            else:
                row["RECLASSIFICATION"] = "PARTIAL"
            writer.writerow(row)


if __name__ == "__main__":
    main()
