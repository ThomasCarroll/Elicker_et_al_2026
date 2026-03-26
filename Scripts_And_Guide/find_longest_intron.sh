#!/bin/bash
set -euo pipefail

# ====================
# User-configurable parameters
# ====================
GENE="R2D2"
FASTA="/ru-auth/local/home/kelicker/SCRATCH/04_Domain_analysis/${GENE}_full_copies.fasta"
GFF_DIR="/ru-auth/local/home/kelicker/SCRATCH/"
OUT="${GENE}_longest_introns.txt"

# ====================
# Start script
# ====================
echo "Processing longest intron of ${GENE}"

declare -A intron_lengths
overall_longest=0
overall_transcript=""

> "$OUT"

# Read all headers
mapfile -t headers < <(grep '^>' "$FASTA")

for header in "${headers[@]}"; do
    header=${header#>}  # remove >
    IFS='_' read -r -a parts <<< "$header"
    n=${#parts[@]}

    transcript="${parts[$((n-2))]}_${parts[$((n-1))]}"
    
    species="${parts[0]}"
    for ((i=1; i<n-2; i++)); do
        species+="_${parts[i]}"
    done

    gff_file="$GFF_DIR/${species}_annotation.gff"
    if [[ ! -f "$gff_file" ]]; then
        echo "GFF file not found for $species" >&2
        continue
    fi

    # Find longest intron
    longest=0
    while read -r len; do
        (( len > longest )) && longest=$len
    done < <(
        awk -F'\t' -v tid="$transcript" '$3=="intron" {
            split($9,a,";");
            parent=""
            for(i in a) if(a[i] ~ /^Parent=/) {split(a[i],b,"="); parent=b[2]}
            if(parent==tid) print $5-$4+1
        }' "$gff_file"
    )

    intron_lengths["$transcript"]=$longest

    # Update overall longest
    if (( longest > overall_longest )); then
        overall_longest=$longest
        overall_transcript=$transcript
    fi

    # Append to output
    echo -e "${transcript}\t${longest}" >> "$OUT"
done

# Append overall longest intron
echo -e "\nOverall longest intron: ${overall_longest} (Transcript: ${overall_transcript})" >> "$OUT"

# ====================
# Calculate Q3 + 3*IQR
# ====================

# Extract intron lengths (skip summary lines)
mapfile -t lengths < <(awk 'NF==2 {print $2}' "$OUT" | sort -n)

quantile() {
    local q=$1
    shift
    local arr=("$@")  # copy all remaining args into array
    local n=${#arr[@]}
    if (( n == 0 )); then
        echo 0
        return
    fi
    local f_idx=$(awk -v n=$n -v q=$q 'BEGIN{print (n-1)*q}')
    local idx_int=${f_idx%.*}
    local frac=$(awk -v f=$f_idx 'BEGIN{print f - int(f)}')
    local val1=${arr[idx_int]}
    local val2=${arr[idx_int+1]:-${arr[idx_int]}}
    awk -v v1=$val1 -v v2=$val2 -v frac=$frac 'BEGIN{printf("%.0f", v1 + (v2-v1)*frac)}'
}

Q1=$(quantile 0.25 "${lengths[@]}")
Q3=$(quantile 0.75 "${lengths[@]}")
IQR=$(( Q3 - Q1 ))
threshold=$(( Q3 + 3*IQR ))

# Append threshold to output
echo -e "\nExtreme outlier threshold (Q3 + 3*IQR): $threshold" >> "$OUT"
echo "Q1 = $Q1, Q3 = $Q3, IQR = $IQR, Q3 + 3*IQR = $threshold"
