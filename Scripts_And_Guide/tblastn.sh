#!/bin/bash
#SBATCH --job-name=AGO2
#SBATCH --cpus-per-task=4

QUERY="AGO2.fasta"
QUERY_BASENAME="$(basename "$QUERY" .fasta)"
GENOME_DIR="/ru-auth/local/home/kelicker/SCRATCH"
OUTDIR="01_tblastn/${QUERY_BASENAME}"
LOG="${OUTDIR}_log.txt"
COMBINED="${OUTDIR}_combined.tsv"

mkdir -p "$OUTDIR" logs
echo "Starting tblastn run at $(date)" > "$LOG"
echo -e "query\tspecies\tchromosome\tqstart\tqend\tsstart\tsend\tlength\tevalue\tbitscore\tpident\tqcovs\tsseq\tstrand" > "$COMBINED"

for GENOME in "$GENOME_DIR"/*_genome.fasta; do
    BASENAME=$(basename "$GENOME" _genome.fasta)
    DB="${GENOME_DIR}/${BASENAME}_genome"
    OUTFILE="${OUTDIR}/${BASENAME}_tblastn.tsv"

    if [ ! -e "${DB}.nin" ]; then
        echo "[$(date)] Making BLAST DB for $BASENAME" | tee -a "$LOG"
        makeblastdb -in "$GENOME" -dbtype nucl -out "$DB" >> "$LOG" 2>&1
    else
        echo "[$(date)] BLAST DB for $BASENAME already exists" | tee -a "$LOG"
    fi

    echo "[$(date)] Running tblastn on $BASENAME" | tee -a "$LOG"
    tblastn \
        -query "$QUERY" \
        -db "$DB" \
        -outfmt "6 qseqid sseqid qstart qend sstart send length evalue bitscore pident qcovs sseq" 2>> "$LOG" | \
    awk -v species="$BASENAME" 'BEGIN { OFS="\t" }
    {
        strand = ($6 < $5) ? "-" : "+"
        print $1, species, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, strand
    }' >> "$COMBINED"

    echo "[$(date)] Finished $BASENAME" | tee -a "$LOG"
done

echo "All tblastn runs complete at $(date)" | tee -a "$LOG"
