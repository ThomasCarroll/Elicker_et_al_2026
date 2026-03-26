#!/usr/bin/env bash
set -euo pipefail

GENE="R2D2"

# === Thresholds per gene ===
declare -A THRESHOLDS=(
  [AGO1]=6289
  [AGO2]=3114
  [AGO3]=6130
  [AUB]=3618
  [DICER1]=257
  [DICER2]=3284
  [DROSHA]=92
  [LOQS]=2972
  [PASHA]=1255
  [PIWI]=7801
  [R2D2]=452
)

# Get threshold for this gene
THRESHOLD="${THRESHOLDS[$GENE]:-}"
if [[ -z "$THRESHOLD" ]]; then
    echo "❌ No threshold defined for gene: $GENE"
    echo "Available genes: ${!THRESHOLDS[@]}"
    exit 1
fi

echo "🔧 Using threshold $THRESHOLD bp for $GENE"

# === Inputs ===
FASTA="/ru-auth/local/home/kelicker/SCRATCH/04_Domain_analysis/${GENE}_partial_copies.fasta"
DIR="/ru-auth/local/home/kelicker/SCRATCH/"

# === Outputs ===
LIKELYFULL="${DIR}/05_Check_contig_ends/${GENE}_partial_copies_annotated_likelyfull.fasta"
REMAINING="${DIR}/05_Check_contig_ends/${GENE}_partial_copies_annotated_remaining.fasta"
SUMMARY="${DIR}/05_Check_contig_ends/${GENE}_partial_copies_annotated_summary.tsv"

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

    # Extract gene ID (last chunk like Aminor_00000013235)
    gene_id=$(echo "$header" | awk -F'_' '{print $(NF-1)"_"$NF}')

    # Get species prefix by removing last two underscore chunks
    species=$(echo "$header" | awk -F'_' '{for (i=1; i<=NF-2; i++) printf (i<NF-2 ? $i"_" : $i)}')

    GFF="${DIR}/${species}_annotation.gff"
    FAI="${DIR}/${species}_genome.fasta.fai"

    # Skip if missing files
    if [[ ! -f "$GFF" || ! -f "$FAI" ]]; then
        echo "⚠️ Missing files for $species — skipping $gene_id"
        continue
    fi

    # Get coordinates
    info=$(get_gene_info "$GFF" "$gene_id") || true
    [[ -z "$info" ]] && continue

    contig=$(echo "$info" | awk '{print $1}')
    start=$(echo "$info" | awk '{print $2}')
    end=$(echo "$info" | awk '{print $3}')
    contig_len=$(get_contig_len "$FAI" "$contig")
    [[ -z "$contig_len" ]] && continue

    near_start=$(( start <= THRESHOLD ? 1 : 0 ))
    near_end=$(( contig_len - end <= THRESHOLD ? 1 : 0 ))

    if (( near_start || near_end )); then
        echo ">$header" >> "$LIKELYFULL"
        echo "$seq" >> "$LIKELYFULL"
        echo -e "${species}\t${gene_id}\t${contig}\t${start}\t${end}\t${contig_len}\tYES" >> "$SUMMARY"
    else
        echo ">$header" >> "$REMAINING"
        echo "$seq" >> "$REMAINING"
        echo -e "${species}\t${gene_id}\t${contig}\t${start}\t${end}\t${contig_len}\tNO" >> "$SUMMARY"
    fi
done

echo "✅ Done!"
echo "Sequences near contig ends → $LIKELYFULL"
echo "Remaining sequences → $REMAINING"
echo "Summary table → $SUMMARY"
