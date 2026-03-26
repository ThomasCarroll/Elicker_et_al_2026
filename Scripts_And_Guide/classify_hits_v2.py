#!/usr/bin/env python3
"""
Parse hmmscan output, classify domain architectures per gene, summarize per species.

FEATURES:
- Flexible domain requirements (either/or, mixed groups)
- True non-overlapping filtering of Ribonuclease domains using coordinates
- Custom rules for DICER1, DICER2, and DROSHA
- Writes a domain architecture table with coordinates

Filters by independent e-value <= 1e-5.
"""

import sys
import csv
from collections import defaultdict, Counter

# -------------------------
# Arguments
# -------------------------
if len(sys.argv) != 5:
    print("Usage: python classify_hits_v2.py <GENE> <domtblout> <classification_out> <species_list>")
    sys.exit(1)

GENE = sys.argv[1].upper()
DOM_TBL = sys.argv[2]
CLASS_OUT = sys.argv[3]
SPECIES_LIST_FILE = sys.argv[4]

IEVALUE_THRESHOLD = 1e-5

ARCH_OUT = f"{GENE}_domain_architecture.tsv"

# -----------------------------------------------------
# Domain definitions
# -----------------------------------------------------

GENE_DOMAINS = {
    "DICER1": {
        "required": {
            "Dicer_PBD": 1,
            "Helicase_C": 1,
            "Dicer_dimer": 1,
            "Dicer_platform": 1,
            "PAZ": 1,
        },
        "mixed_groups": [
            {
                "domains": ["Ribonuclease_3", "Ribonuclease_3_3", "Ribonucleas_3_3"],
                "count": 2
            }
        ]
    },

    "DICER2": {
        "required": {
            "Dicer_dimer": 1,
            "Dicer_platform": 1,
            "PAZ": 1,
            "Helicase_C": 1,
            "Dicer_dsRBD": 1
        },
        "either_or": [
            ["DEAD", "ResIII"]
        ],
        "mixed_groups": [
            {
                "domains": ["Ribonuclease_3", "Ribonuclease_3_3", "Ribonucleas_3_3"],
                "count": 2
            }
        ]
    },

    "DROSHA": {
        "required": {
            "dsrm": 1
        },
        "mixed_groups": [
            {
                "domains": ["Ribonuclease_3", "Ribonuclease_3_3", "Ribonucleas_3_3"],
                "count": 2
            }
        ]
    },

    "PASHA": {"required": {"dsrm": 1}},
    "LOQS": {"required": {"dsrm": 2}},
    "R2D2": {"required": {"dsrm": 2}},
    "AGO1": {"required": {"Piwi": 1, "ArgoMid": 1, "ArgoN": 1, "PAZ": 1, "ArgoL1": 1, "ArgoL2": 1}},
    "AGO2": {"required": {"Piwi": 1, "ArgoMid": 1, "ArgoN": 1, "PAZ": 1, "ArgoL1": 1, "ArgoL2": 1}},
    "AGO3": {"required": {"Piwi": 1, "PAZ": 1, "Piwi_N": 1}},
    "PIWI": {"required": {"Piwi": 1, "PAZ": 1, "Piwi_N": 1}},
    "AUB": {"required": {"Piwi": 1, "ArgoN": 1, "PAZ": 1}}
}

if GENE not in GENE_DOMAINS:
    print(f"Error: Unknown gene '{GENE}'")
    sys.exit(1)

EXPECTED = GENE_DOMAINS[GENE]
REQUIRED = EXPECTED.get("required", {})
EITHER_OR = EXPECTED.get("either_or", [])
MIXED = EXPECTED.get("mixed_groups", [])

# -------------------------
# Load species list
# -------------------------
ALL_SPECIES = []
with open(SPECIES_LIST_FILE) as f:
    for line in f:
        sp = line.strip()
        if sp:
            ALL_SPECIES.append(sp)

# ------------------------------------------------------
# Parse hmmscan domtblout with coordinates
# ------------------------------------------------------

dom_hits = defaultdict(list)

with open(DOM_TBL) as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.split()
        if len(fields) < 23:
            continue

        target_name = fields[0]
        query_name = fields[3]

        try:
            i_evalue = float(fields[12])
            ali_from = int(fields[17])
            ali_to = int(fields[18])
        except ValueError:
            continue

        if i_evalue <= IEVALUE_THRESHOLD:
            dom_hits[query_name].append({
                "domain": target_name,
                "start": ali_from,
                "end": ali_to
            })

# ------------------------------------------------------
# Helper function: remove overlapping domains
# ------------------------------------------------------

def filter_non_overlapping(hits, allowed_domains):
    relevant = [h for h in hits if h["domain"] in allowed_domains]
    relevant.sort(key=lambda x: (x["end"] - x["start"]), reverse=True)

    kept = []
    for h in relevant:
        if all(h["end"] < k["start"] or h["start"] > k["end"] for k in kept):
            kept.append(h)
    return kept

# ------------------------------------------------------
# Classification
# ------------------------------------------------------

RNASE_DOMAINS = ["Ribonuclease_3", "Ribonuclease_3_3", "Ribonucleas_3_3"]

with open(CLASS_OUT, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["query", "species", "gene", "domains_detected", "domain_counts", "classification"])

    for q in sorted(dom_hits):
        hits = dom_hits[q]

        # RNase domains: non-overlapping
        non_overlapping = filter_non_overlapping(hits, RNASE_DOMAINS)

        domain_counts = Counter()
        for h in hits:
            if h["domain"] not in RNASE_DOMAINS:
                domain_counts[h["domain"]] += 1
        for h in non_overlapping:
            domain_counts[h["domain"]] += 1

        ### FIX: detection is independent of rule satisfaction
        found_any = len(domain_counts) > 0

        ### FIX: rule checks ONLY affect is_full
        is_full = True

        for dom, required_count in REQUIRED.items():
            if domain_counts.get(dom, 0) < required_count:
                is_full = False

        for options in EITHER_OR:
            if not any(domain_counts.get(d, 0) > 0 for d in options):
                is_full = False

        for group in MIXED:
            total = sum(domain_counts.get(d, 0) for d in group["domains"])
            if total < group["count"]:
                is_full = False

        if is_full:
            classification = "FULL"
        elif found_any:
            classification = "PARTIAL"
        else:
            classification = "IGNORE"

        parts = q.split("_")
        species = "_".join(parts[:-2])

        count_str = ",".join(f"{d}({domain_counts[d]})" for d in sorted(domain_counts))

        w.writerow([
            q,
            species,
            GENE,
            ",".join(sorted(domain_counts)),
            count_str,
            classification
        ])

# ------------------------------------------------------
# Species summary
# ------------------------------------------------------

summary = {sp: {"FULL": 0, "PARTIAL": 0, "IGNORE": 0} for sp in ALL_SPECIES}

with open(CLASS_OUT) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if row["species"] in summary:
            summary[row["species"]][row["classification"]] += 1

with open(f"{GENE}_species_summary.tsv", "w", newline="") as out:
    w = csv.writer(out, delimiter="\t")
    w.writerow(["species", "FULL", "PARTIAL", "IGNORE"])
    for sp in ALL_SPECIES:
        w.writerow([sp, summary[sp]["FULL"], summary[sp]["PARTIAL"], summary[sp]["IGNORE"]])
