#!/bin/bash
#SBATCH --job-name=PIWI_reciprocal
#SBATCH --cpus-per-task=4

set -o pipefail

###############################
#       CHANGE THIS ONLY      #
GENE="PIWI"
GENE_ID="piwi"
###############################

# Input files
QUERY_FASTA="${GENE}_predicted_proteins_clean.fa"
BLAST_DB="Drosophila_melanogaster_genome"
GFF_FILE="Drosophila_melanogaster_annotation.gff"

# Output files
CONFIRMED="${GENE}_NOCDS_AUGUSTUS_confirmed.fasta"
ELIMINATED="${GENE}_NOCDS_AUGUSTUS_eliminated.fasta"
TMP_DIR="${GENE}_reciprocal_blast_tmp"

# Setup
mkdir -p "$TMP_DIR"
: > "$CONFIRMED"
: > "$ELIMINATED"

# Counters
total_queries=0
confirmed_count=0
eliminated_count=0

# Build header–sequence table
awk '/^>/ {if(seq) print seq; print; seq=""; next} {seq = seq $0} END {print seq}' "$QUERY_FASTA" \
    | paste - - > "$TMP_DIR/one_line_queries.tsv"

while IFS=$'\t' read -r header sequence; do
    ((total_queries++))
    id=$(echo "$header" | sed 's/^>//;s/ .*//')
    tmp_query="$TMP_DIR/${id}.fasta"
    printf "%s\n%s\n" "$header" "$sequence" > "$tmp_query"

    # Run BLAST
    blast_output="$TMP_DIR/${id}_blast.txt"
    tblastn -query "$tmp_query" -db "$BLAST_DB" -evalue 1e-5 \
            -outfmt "6 sseqid sstart send" -max_target_seqs 1 \
            > "$blast_output" 2>/dev/null

    # No hit
    if [ ! -s "$blast_output" ]; then
        echo "No hit for $id"
        printf "%s\n%s\n" "$header" "$sequence" >> "$ELIMINATED"
        ((eliminated_count++))
        continue
    fi

    read -r contig start end < <(head -n1 "$blast_output")

    # Ensure start <= end
    if [ "$start" -gt "$end" ]; then
        tmp=$start; start=$end; end=$tmp
    fi

    # Find overlapping gene
    gene_name=$(awk -v contig="$contig" -v start="$start" -v end="$end" '
        $1 == contig && $3 == "gene" && $4 <= end && $5 >= start {
            if (match($0, /[Gg]ene=([^;]+)/, arr))
                print arr[1]
        }
    ' "$GFF_FILE" | head -n 1)

    # Classification — automatically matches the gene name
    if [[ -n "$gene_name" && "$gene_name" == *"$GENE_ID"* ]]; then
        echo "$id --> $gene_name (Confirmed)"
        printf "%s\n%s\n" "$header" "$sequence" >> "$CONFIRMED"
        ((confirmed_count++))
    else
        echo "$id --> ${gene_name:-<no_gene>} (Eliminated)"
        printf "%s\n%s\n" "$header" "$sequence" >> "$ELIMINATED"
        ((eliminated_count++))
    fi
done < "$TMP_DIR/one_line_queries.tsv"

# Report
echo
echo "===== Reciprocal BLAST Summary for $GENE ====="
echo "Total queries:    $total_queries"
echo "Confirmed:        $confirmed_count"
echo "Eliminated:       $eliminated_count"
