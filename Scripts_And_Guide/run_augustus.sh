#!/bin/bash
# Gene prediction pipeline with reciprocal BLAST, HMMER, and classification
# Usage: ./gene_pipeline.sh GENE_NAME [REFERENCE_GENE]
# Example: ./gene_pipeline.sh DICER Dcr-1

#SBATCH --cpus-per-task=4
set -euo pipefail

# ------------------------------
# Input gene name
# ------------------------------
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 GENE_NAME"
  exit 1
fi

GENE="$1"

echo "Running Augustus prediction for gene: $GENE"

# ------------------------------
# Input & Parameters
# ------------------------------
BED_FILE="${GENE}_no_missingGFF.bed"
EXTEND_BP=500
GENOME_DIR="../"
OUTPUT_DIR="${GENE}_gene_prediction"
PFAM_DB="../04_Domain_analysis/Pfam-A.hmm"
AUGUSTUS_SPECIES="fly"

BLAST_DB="/ru-auth/local/home/kelicker/SCRATCH/Drosophila_melanogaster_genome"
GFF_FILE="/ru-auth/local/home/kelicker/SCRATCH/Drosophila_melanogaster_annotation.gff"

mkdir -p "$OUTPUT_DIR"
TMP_DIR="$OUTPUT_DIR/tmp"
mkdir -p "$TMP_DIR"

command -v bedtools >/dev/null
command -v samtools >/dev/null
command -v augustus >/dev/null
command -v tblastn >/dev/null
command -v hmmscan >/dev/null

# ------------------------------
# 0. Fix BED - skip header and include strand
# ------------------------------
awk 'BEGIN{OFS="\t"} NR>1 && NF>=7 {s=$2; e=$3; if(s>e){t=s; s=e; e=t} print $1,s,e,$4,$5,$7}' \
  "$BED_FILE" > "$OUTPUT_DIR/${GENE}_fixed.bed"

# ------------------------------
# 1. Sort BED
# ------------------------------
bedtools sort -i "$OUTPUT_DIR/${GENE}_fixed.bed" > "$OUTPUT_DIR/${GENE}_sorted.bed"

# ------------------------------
# 2. Extend regions
# ------------------------------
: > "$OUTPUT_DIR/${GENE}_extended.bed"
while read -r chrom start end species score strand; do
  genome_fasta="${GENOME_DIR}/${species}_genome.fasta"
  [[ ! -f $genome_fasta ]] && echo "WARNING: $genome_fasta not found, skipping" && continue
  [[ ! -f ${genome_fasta}.fai ]] && samtools faidx "$genome_fasta"
  contig_len=$(awk -v c="$chrom" '$1==c{print $2}' "${genome_fasta}.fai")
  [[ -z "${contig_len:-}" ]] && echo "WARNING: Contig $chrom not found" && continue
  new_start=$(( start - EXTEND_BP )); (( new_start < 0 )) && new_start=0
  new_end=$(( end + EXTEND_BP )); (( new_end > contig_len )) && new_end=$contig_len
  printf "%s\t%d\t%d\t%s\t%s\n" "$chrom" "$new_start" "$new_end" "$species" "$strand" \
    >> "$OUTPUT_DIR/${GENE}_extended.bed"
done < "$OUTPUT_DIR/${GENE}_sorted.bed"

# ------------------------------
# 3. Extract sequences
# ------------------------------
: > "$OUTPUT_DIR/${GENE}_sequences.fa"
while read -r chrom start end species strand; do
  genome_fasta="${GENOME_DIR}/${species}_genome.fasta"
  [[ ! -f $genome_fasta ]] && continue
  tmp_bed=$(mktemp)
  printf "%s\t%d\t%d\t%s|%s:%d-%d(%s)\t0\t%s\n" "$chrom" "$start" "$end" "$species" "$chrom" "$start" "$end" "$strand" "$strand" > "$tmp_bed"
  bedtools getfasta -fi "$genome_fasta" -bed "$tmp_bed" -s -name >> "$OUTPUT_DIR/${GENE}_sequences.fa"
  rm -f "$tmp_bed"
done < "$OUTPUT_DIR/${GENE}_extended.bed"
echo "Sequences extracted: $OUTPUT_DIR/${GENE}_sequences.fa"

# ------------------------------
# 4. Run AUGUSTUS
# ------------------------------
AUG_OUT="$OUTPUT_DIR/${GENE}_augustus_output.txt"
augustus --species="$AUGUSTUS_SPECIES" "$OUTPUT_DIR/${GENE}_sequences.fa" --gff3=on --protein=on > "$AUG_OUT"
echo "AUGUSTUS prediction complete: $AUG_OUT"

# ------------------------------
# 5. Clean AUGUSTUS output and extract protein sequences
# ------------------------------

CLEAN_PROT="$OUTPUT_DIR/${GENE}_predicted_proteins_clean.fa"

awk '
BEGIN { seq=""; header=""; inseq=0 }

/^# ----- prediction on sequence number/ {
    if (seq != "" && header != "") {
        # remove trailing "Evidence for and against this transcript:" if present
        sub(/Evidence for and against this transcript:.*$/,"", seq)
        print ">"header "\n" seq
        seq=""; header=""; inseq=0
    }
}

/name =/ {
    if (match($0,/name = ([^:]+:[0-9]+-[0-9]+\([+-]\))/ ,a)) {
        header=a[1]
    }
}

/^# protein sequence = \[/ {
    inseq=1
    seq=""
    sub(/^# protein sequence = \[/,"")
    sub(/\].*$/,"")
    seq=$0
    next
}

/^#/ && inseq {
    gsub(/^# */,"")
    if ($0 ~ /\]$/ || $0 ~ /^Evidence for and against/) {
        sub(/\].*$/,"")
        seq=seq $0
        # remove trailing "Evidence for and against this transcript:" if it somehow slipped in
        sub(/Evidence for and against this transcript:.*$/,"", seq)
        inseq=0
        if(header!="") print ">"header "\n" seq
        seq=""; header=""
    } else {
        seq=seq $0
    }
}

END {
    if(seq != "" && header != "") {
        sub(/Evidence for and against this transcript:.*$/,"", seq)
        print ">"header "\n" seq
    }
}
' "$AUG_OUT" > "$CLEAN_PROT"

echo "Cleaned protein sequences: $CLEAN_PROT"