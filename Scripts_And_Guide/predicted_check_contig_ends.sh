#!/usr/bin/env bash
set -euo pipefail

GENE="R2D2"

# === Inputs ===
FASTA="${GENE}_NOCDS_AUGUSTUS_confirmed_partial_copies.fasta"
DIR="/ru-auth/local/home/kelicker/SCRATCH/"

# ------------------------------
# Contig-end thresholds
# ------------------------------
declare -A THRESHOLDS=( 
    ["AGO1"]=6289 ["AGO2"]=3114 ["AGO3"]=6130 ["AUB"]=3618 
    ["DICER1"]=257 ["DICER2"]=3284 ["DROSHA"]=92 ["LOQS"]=2972 
    ["PASHA"]=1255 ["PIWI"]=7801 ["R2D2"]=452 
)

# Use the appropriate threshold for this gene
THRESHOLD=${THRESHOLDS[$GENE]:-500}  # default 500 if not defined

# === Outputs ===
LIKELYFULL="${GENE}_NOCDS_AGUSTUS_confirmed_partial_likelyfull.fasta"
REMAINING="${GENE}_NOCDS_AGUSTUS_confirmed_partial_remaining.fasta"
SUMMARY="${GENE}_NOCDS_AGUSTUS_confirmed_partial_summary.tsv"

> "$LIKELYFULL"
> "$REMAINING"
echo -e "species\tgene_id\tcontig\tstart\tend\tcontig_len\tnear_end" > "$SUMMARY"

# === Helper functions ===
get_gene_info() {
    local gff="$1"
    local gene_id="$2"
    awk -v id="$gene_id" -F '\t' '$3=="mRNA" && $9~id {print $1, $4, $5; exit}' "$gff"
}

get_contig_len() {
    local fai="$1"
    local contig="$2"
    awk -v c="$contig" '$1==c {print $2; exit}' "$fai"
}

# === Process each FASTA entry ===
awk '/^>/ {if (seq) print header"\t"seq; header=$0; seq=""} /^[^>]/ {seq=seq$0} END {print header"\t"seq}' "$FASTA" |
while IFS=$'\t' read -r header seq; do
    # Remove leading ">" if present
    header=${header#>}

    # Extract species (everything before the first |)
    species="${header%%|*}"

    # Extract contig and coordinates from header (after |)
    contig_info="${header#*|}"  # CM009976.1:5167072-5168524(+)
    if [[ $contig_info =~ ^([^:]+):([0-9]+)-([0-9]+)\(([+-])\)$ ]]; then
        contig="${BASH_REMATCH[1]}"
        start="${BASH_REMATCH[2]}"
        end="${BASH_REMATCH[3]}"
        strand="${BASH_REMATCH[4]}"
    else
        echo "⚠️ Header format unexpected: $header — skipping"
        continue
    fi

    FAI="${DIR}/${species}_genome.fasta.fai"

    # Skip if missing file
    if [[ ! -f "$FAI" ]]; then
        echo "⚠️ Missing FAI for $species — skipping $header"
        continue
    fi

    contig_len=$(get_contig_len "$FAI" "$contig")
    [[ -z "$contig_len" ]] && continue

    # --- Use appropriate threshold per gene ---
    THRESHOLD=${THRESHOLDS[$GENE]:-500}

    near_start=$(( start <= THRESHOLD ? 1 : 0 ))
    near_end=$(( contig_len - end <= THRESHOLD ? 1 : 0 ))

    if (( near_start || near_end )); then
        echo ">$header" >> "$LIKELYFULL"
        echo "$seq" >> "$LIKELYFULL"
        echo -e "${species}\t${header}\t${contig}\t${start}\t${end}\t${contig_len}\tYES" >> "$SUMMARY"
    else
        echo ">$header" >> "$REMAINING"
        echo "$seq" >> "$REMAINING"
        echo -e "${species}\t${header}\t${contig}\t${start}\t${end}\t${contig_len}\tNO" >> "$SUMMARY"
    fi
done

echo "✅ Done!"
echo "Sequences near contig ends → $LIKELYFULL"
echo "Remaining sequences → $REMAINING"
echo "Summary table → $SUMMARY"
