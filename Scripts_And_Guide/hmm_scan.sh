#!/bin/bash
#SBATCH --cpus-per-task=8

# -------------------------
# Usage: sbatch hmm_scan.sh AGO1
# -------------------------
if [ $# -ne 1 ]; then
    echo "Usage: sbatch $0 <GENE>"
    exit 1
fi

GENE=$1

CANDIDATES="${GENE}_newhits_filtered_confirmed.fasta"
PFAM_DB="Pfam-A.hmm"
DOM_TBL="${GENE}_domtblout.tsv"
CLASS_OUT="${GENE}_classification.tsv"
SPECIES_LIST="species_list.txt"
SPECIES_SUM="${GENE}_species_summary.tsv"

FULL_FILE="${GENE}_full_copies.fasta"
PARTIAL_FILE="${GENE}_partial_copies.fasta"
IGNORE_FILE="${GENE}_ignore_copies.fasta"

echo "======================================"
echo "Running hmmscan for $GENE..."
echo "======================================"

hmmscan --cpu 8 --domtblout "$DOM_TBL" "$PFAM_DB" "$CANDIDATES" > "${GENE}_hmmscan.log"

echo "======================================"
echo "Classifying domains for $GENE..."
echo "======================================"

python classify_hits_v2.py "$GENE" "$DOM_TBL" "$CLASS_OUT" "$SPECIES_LIST"

echo "======================================"
echo "Extracting FASTA sequences by class..."
echo "======================================"

# Run the output FASTA logic *inside* this script
> "$FULL_FILE"
> "$PARTIAL_FILE"

# Counters
full_count=0
partial_count=0
ignore_count=0
skipped_count=0

# Load classification table
declare -A class_map
while IFS=$'\t' read -r query species gene domains domain_counts classification; do
    [[ $query == "query" ]] && continue
    classification=$(echo "$classification" | tr -d '\r' | xargs)
    class_map["$query"]="$classification"
done < "$CLASS_OUT"

# Parse FASTA and categorize
seq_name=""
seq=""
while read -r line; do
    if [[ $line == ">"* ]]; then
        if [[ -n $seq_name ]]; then
            class=${class_map[$seq_name]}
            if [[ $class == "FULL" ]]; then
                echo -e ">$seq_name\n$seq" >> "$FULL_FILE"
                ((full_count++))
            elif [[ $class == "PARTIAL" ]]; then
                echo -e ">$seq_name\n$seq" >> "$PARTIAL_FILE"
                ((partial_count++))
            elif [[ $class == "IGNORE" ]]; then
                echo -e ">$seq_name\n$seq" >> "$IGNORE_FILE"
                ((ignore_count++))
            else
                ((skipped_count++))
            fi
        fi
        seq_name="${line#>}"
        seq=""
    else
        seq+="$line"$'\n'
    fi
done < "$CANDIDATES"

# Last sequence
if [[ -n $seq_name ]]; then
    class=${class_map[$seq_name]}
    if [[ $class == "FULL" ]]; then
        echo -e ">$seq_name\n$seq" >> "$FULL_FILE"
        ((full_count++))
    elif [[ $class == "PARTIAL" ]]; then
        echo -e ">$seq_name\n$seq" >> "$PARTIAL_FILE"
        ((partial_count++))
    else
        ((skipped_count++))
    fi
fi

echo "======================================"
echo "Pipeline complete for $GENE!"
echo "======================================"
echo "Files created:"
echo "  - domtblout:          $DOM_TBL"
echo "  - classification:     $CLASS_OUT"
echo "  - species summary:    $SPECIES_SUM"
echo "  - FULL sequences:     $FULL_FILE"
echo "  - PARTIAL sequences:  $PARTIAL_FILE"
echo "  - IGNORE sequences:   $IGNORE_FILE"
echo
echo "Counts:"
echo "  FULL:      $full_count"
echo "  PARTIAL:   $partial_count"
echo "  IGNORE:    $ignore_count"
echo "  SKIPPED:   $skipped_count"
