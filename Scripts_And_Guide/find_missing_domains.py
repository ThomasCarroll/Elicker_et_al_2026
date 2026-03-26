#!/usr/bin/env python3

import sys, os, subprocess
from collections import Counter

if len(sys.argv) != 2:
    print("Usage: python classify_missing_domains.py <GENE>")
    sys.exit(1)

GENE = sys.argv[1].upper()
SCRATCH = "/ru-auth/local/home/kelicker/SCRATCH"

INPUT_FILE = f"{SCRATCH}/07_Find_missing_domains/{GENE}_record_keeping_with_domains_with_contig_end.txt"
OUTPUT_FILE = f"{SCRATCH}/07_Find_missing_domains/{GENE}_record_keeping_missing_domains.txt"
PFAM = f"{SCRATCH}/04_Domain_analysis/Pfam-A.hmm"

EXTEND = 5000
os.makedirs(f"{GENE}_hmmer", exist_ok=True)

# -------------------------
# Domain architecture rules
# -------------------------
REQUIRED = {}
EITHER_OR = []
MIXED = []

if GENE in ["AGO1","AGO2"]:
    REQUIRED = {"Piwi":1,"ArgoMid":1,"ArgoN":1,"PAZ":1,"ArgoL1":1,"ArgoL2":1}
elif GENE in ["AGO3","PIWI"]:
    REQUIRED = {"Piwi":1,"PAZ":1,"Piwi_N":1}
elif GENE=="AUB":
    REQUIRED = {"Piwi":1,"ArgoN":1,"PAZ":1}
elif GENE=="PASHA":
    REQUIRED = {"dsrm":1}
elif GENE in ["LOQS","R2D2"]:
    REQUIRED = {"dsrm":2}
elif GENE=="DROSHA":
    REQUIRED={"dsrm":1}
    MIXED=[{"domains":["Ribonuclease_3","Ribonuclease_3_3","Ribonucleas_3_3"],"count":2}]
elif GENE=="DICER1":
    REQUIRED={"Dicer_PBD":1,"Helicase_C":1,"Dicer_dimer":1,"Dicer_platform":1,"PAZ":1,"Dicer_dsRBD":1}
    MIXED=[{"domains":["Ribonuclease_3","Ribonuclease_3_3","Ribonucleas_3_3"],"count":2}]
elif GENE=="DICER2":
    REQUIRED={"Dicer_dimer":1,"Dicer_platform":1,"PAZ":1,"Helicase_C":1,"Dicer_dsRBD":1}
    EITHER_OR=[["DEAD","ResIII"]]
    MIXED=[{"domains":["Ribonuclease_3","Ribonuclease_3_3","Ribonucleas_3_3"],"count":2}]
else:
    print(f"Unknown gene: {GENE}")
    sys.exit(1)

# -------------------------
# Load input table
# -------------------------
hits = []
with open(INPUT_FILE) as f:
    header = f.readline().strip().split("\t")
    for line in f:
        if not line.strip():
            continue
        hits.append(dict(zip(header,line.strip().split("\t"))))

# -------------------------
# Process hits
# -------------------------
with open(OUTPUT_FILE,"w") as out:
    out.write("\t".join(header + ["FOUND_DOMAINS","RECLASSIFICATION"]) + "\n")

    for row in hits:

        original_class = row.get("CLASSIFICATION",".")

        # Preserve FULL genes
        if original_class == "FULL":
            out.write("\t".join([row[h] for h in header]+[".","FULL"])+"\n")
            continue

        # Pass through non-PARTIAL
        if original_class != "PARTIAL":
            out.write("\t".join([row[h] for h in header]+[".",original_class])+"\n")
            continue

        species = row["SPECIES"]
        contig = row["CHROMOSOME"]
        start = int(row["START"])
        end = int(row["END"])

        genome_file = f"{SCRATCH}/{species}_genome.fasta"
        fai_file = genome_file + ".fai"

        # Get contig length
        contig_len = None
        with open(fai_file) as f:
            for l in f:
                parts = l.strip().split("\t")
                if parts[0] == contig:
                    contig_len = int(parts[1])
                    break

        if contig_len is None:
            out.write("\t".join([row[h] for h in header]+[".","UNKNOWN_CONTIG"])+"\n")
            continue

        ext_start = max(1,start-EXTEND)
        ext_end = min(contig_len,end+EXTEND)

        region_fasta = f"{GENE}_hmmer/{species}_{row['HIT_ID']}_region.fa"
        region_6fa = f"{GENE}_hmmer/{species}_{row['HIT_ID']}_6frame.fa"
        domtblout = f"{GENE}_hmmer/{species}_{row['HIT_ID']}_domtbl.txt"

        # Extract region
        subprocess.run(
            f"samtools faidx {genome_file} {contig}:{ext_start}-{ext_end} > {region_fasta}",
            shell=True, check=True
        )

        # 6-frame translation
        subprocess.run(
            ["transeq","-sequence",region_fasta,"-outseq",region_6fa,"-frame","6"],
            check=True
        )

        # hmmscan
        subprocess.run(
            ["hmmscan","--domtblout",domtblout,PFAM,region_6fa],
            check=True
        )

        # -------------------------
        # Parse HMMer output
        # -------------------------
        hmmer_counts = Counter()

        with open(domtblout) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.split()
                if len(fields) < 23:
                    continue
                try:
                    full_eval = float(fields[6])
                    dom_ie = float(fields[11])
                except:
                    continue
                if full_eval < 1e-5 and dom_ie < 1e-5:
                    hmmer_counts[fields[0]] += 1

        # -------------------------
        # Parse annotated DOMAINS column
        # -------------------------
        annotated_counts = Counter()

        annotated_field = row.get("DOMAINS","")
        if annotated_field and annotated_field != ".":
            for d in annotated_field.split(","):
                d = d.strip()
                if "(" in d:
                    name = d.split("(")[0]
                    try:
                        count = int(d.split("(")[1].rstrip(")"))
                    except:
                        count = 1
                else:
                    name = d
                    count = 1
                annotated_counts[name] += count

        # -------------------------
        # Architecture relevant domains
        # -------------------------
        arch_domains = set(REQUIRED.keys())

        for group in MIXED:
            arch_domains.update(group["domains"])

        for group in EITHER_OR:
            arch_domains.update(group)

        hmmer_relevant = Counter({
            d: hmmer_counts[d] for d in hmmer_counts if d in arch_domains
        })

        # -------------------------
        # Evaluate architecture rescue
        # -------------------------
        contributed = False
        architecture_satisfied = True

        # REQUIRED
        for dom, req in REQUIRED.items():
            annotated = annotated_counts.get(dom,0)
            hmmer = hmmer_relevant.get(dom,0)

            already = min(req, annotated)
            now_total = min(req, annotated + hmmer)

            if now_total > already:
                contributed = True

            if annotated + hmmer < req:
                architecture_satisfied = False

        # MIXED
        for group in MIXED:
            required = group["count"]
            annotated_total = sum(annotated_counts.get(d,0) for d in group["domains"])
            hmmer_total = sum(hmmer_relevant.get(d,0) for d in group["domains"])

            already = min(required, annotated_total)
            now_total = min(required, annotated_total + hmmer_total)

            if now_total > already:
                contributed = True

            if annotated_total + hmmer_total < required:
                architecture_satisfied = False

        # EITHER/OR
        for group in EITHER_OR:
            annotated_present = any(annotated_counts.get(d,0) > 0 for d in group)
            hmmer_present = any(hmmer_relevant.get(d,0) > 0 for d in group)

            if not annotated_present and hmmer_present:
                contributed = True

            if not (annotated_present or hmmer_present):
                architecture_satisfied = False

        # -------------------------
        # Reclassification
        # -------------------------
        if architecture_satisfied:
            reclass = "LIKELY_FULL_DOMAINS_FOUND"
        elif contributed:
            reclass = "PARTIAL_DOMAINS_FOUND"
        else:
            reclass = "PARTIAL"

        found_display = [
            f"{d}({hmmer_relevant[d]})"
            for d in sorted(hmmer_relevant)
        ]

        out.write("\t".join(
            [row[h] for h in header] +
            [",".join(found_display), reclass]
        )+"\n")

print(f"\nDone. Output written to {OUTPUT_FILE}")