#!/bin/bash

GENE_NAME="R2D2"

BED_FILE="${GENE_NAME}_newhits_filtered.txt"
OUTPUT_FASTA="${GENE_NAME}_newhits_filtered.fasta"
NO_CDS_BED="${GENE_NAME}_newhits_filtered_noCDS.bed"

# Create or clear the output files
> "$OUTPUT_FASTA"
echo -e "contig\tstart\tend\tspecies\thit_id\treason\tstrand" > "$NO_CDS_BED"

echo "Started extracting protein sequences for gene: $GENE_NAME"
echo "Input BED file: $BED_FILE."

# Declare an associative array to track processed mRNAs
declare -A mrna_seen

# Counters
fasta_count=0
no_cds_count=0

# Read BED file line by line, skipping the header
tail -n +2 "$BED_FILE" | while IFS=$'\t' read -r query species hit_id chromosome qstart qend sstart send length evalue bitscore pident qcovs sseq strand hit_percent_identity hit_coverage_percent; do

    contig="$chromosome"
    start="$sstart"
    end="$send"

    echo "Processing region: $contig:$start-$end for species: $species"

    GFF_FILE="/ru-auth/local/home/kelicker/SCRATCH/${species}_annotation.gff"
    GENOME_FASTA="/ru-auth/local/home/kelicker/SCRATCH/${species}_genome.fasta"
    
    echo "Using GFF file: $GFF_FILE"
    echo "Using genome FASTA: $GENOME_FASTA"

    # Check if the GFF file exists
    if [[ ! -f "$GFF_FILE" ]]; then
        echo "Error: GFF file $GFF_FILE not found!" >&2
        echo -e "$contig\t$start\t$end\t$species\t$hit_id\tMissing_GFF\t$strand" >> "$NO_CDS_BED"
        ((no_cds_count++))
        continue
    fi

    # Get all mRNA IDs overlapping the region
    mapfile -t mrna_ids < <(awk -v chr="$contig" -v start="$start" -v end="$end" '
        $1 == chr && $4 <= end && $5 >= start && $3 == "mRNA" {
            match($9, /ID=([^;]+)/, a)
            if (a[1] != "") print a[1]
        }' "$GFF_FILE")
    
    # Choose longest mRNA
    longest_mrna=""
    max_cds_len=0
    
    for m_id in "${mrna_ids[@]}"; do
        cds_len=$(awk -v id="$m_id" '$3 == "CDS" && $9 ~ ("Parent=" id) {sum += $5 - $4 + 1} END {print sum}' "$GFF_FILE")
        if [[ -n "$cds_len" && "$cds_len" -gt "$max_cds_len" ]]; then
            max_cds_len="$cds_len"
            longest_mrna="$m_id"
        fi
    done
    
    mrna_id="$longest_mrna"

    # No mRNA found
    if [[ -z "$mrna_id" ]]; then
        echo -e "$contig\t$start\t$end\t$species\t$hit_id\tNo_mRNA_found\t$strand" >> "$NO_CDS_BED"
        ((no_cds_count++))
        continue
    fi

    full_mrna_id="${species}_${mrna_id}"
    echo "Full mRNA ID: $full_mrna_id"

    # Skip duplicates
    if [[ -n "${mrna_seen[$full_mrna_id]}" ]]; then
        echo "mRNA $full_mrna_id already processed, skipping."
        continue
    fi
    mrna_seen[$full_mrna_id]=1  

    CDS_GFF="temp_${species}_${mrna_id}.gff"
    CDS_FASTA="temp_${species}_${mrna_id}.fasta"

    # Extract CDS entries
    awk -v id="$mrna_id" '$3 == "CDS" && $9 ~ ("Parent=" id) {print}' "$GFF_FILE" > "$CDS_GFF"

    # No CDS found
    if [[ ! -s "$CDS_GFF" ]]; then
        echo -e "$contig\t$start\t$end\t$species\t$hit_id\tNo_CDS_found\t$strand" >> "$NO_CDS_BED"
        ((no_cds_count++))
        rm -f "$CDS_GFF"
        continue
    fi

    # Extract the CDS sequence
    gffread -g "$GENOME_FASTA" -x "$CDS_FASTA" "$CDS_GFF"
    rm -f "$CDS_GFF"

    if [[ ! -s "$CDS_FASTA" ]]; then
        echo -e "$contig\t$start\t$end\t$species\t$hit_id\tCDS_extraction_failed\t$strand" >> "$NO_CDS_BED"
        ((no_cds_count++))
        rm -f "$CDS_FASTA"
        continue
    fi

    # Translate sequence
    PROTEIN_FASTA="temp_${species}_${mrna_id}_protein.fasta"
    transeq -sequence "$CDS_FASTA" -outseq "$PROTEIN_FASTA" -frame=1 -clean
    
    if [[ ! -s "$PROTEIN_FASTA" ]]; then
        echo -e "$contig\t$start\t$end\t$species\t$hit_id\tTranslation_failed\t$strand" >> "$NO_CDS_BED"
        ((no_cds_count++))
        rm -f "$CDS_FASTA" "$PROTEIN_FASTA"
        continue
    fi
    
    # Add protein to output FASTA
    echo ">$species"_"$mrna_id" >> "$OUTPUT_FASTA"
    sed '1d' "$PROTEIN_FASTA" >> "$OUTPUT_FASTA"
    ((fasta_count++))
    
    rm -f "$CDS_FASTA" "$PROTEIN_FASTA"

done

echo "Processing complete!"
echo "$fasta_count annotated sequences saved in FASTA format at $OUTPUT_FASTA."
echo "$no_cds_count regions without CDS or protein saved in $NO_CDS_BED."
