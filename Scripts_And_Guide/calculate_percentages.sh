#!/bin/bash
set -euo pipefail

# ------------------------------
# Input file and basename
# ------------------------------
RESULT="AGO2_combined.tsv"
BASENAME=$(basename "$RESULT" _combined.tsv)
echo "Processing $RESULT"

# ------------------------------
# Sort by species, chromosome, sstart
# ------------------------------
SORTED="${BASENAME}_combined.tsv.sorted"
GROUPED="${BASENAME}_combined.tsv.sorted.grouped"

echo "Sorting $RESULT..."
head -n 1 "$RESULT" > "$SORTED"
tail -n +2 "$RESULT" | sort -t$'\t' -k2,2 -k3,3 -k6,6n >> "$SORTED"

# ------------------------------
# Assign hit_id
# ------------------------------
echo "Assigning hit_id..."
awk -F'\t' -v OFS='\t' '
BEGIN { hit_id = 0 }
NR==1 { print $1, $2, "hit_id", $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14; next }
{
    species=$2; chrom=$3; strand=$14; qstart=$4+0

    if (species != prev_species || chrom != prev_chrom || strand != prev_strand) { hit_id++ }
    else {
        if(strand=="+"){ if(qstart<=prev_qstart) hit_id++ }
        else if(strand=="-"){ if(qstart>=prev_qstart) hit_id++ }
    }

    print $1, $2, hit_id, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14

    prev_species=species; prev_chrom=chrom; prev_strand=strand; prev_qstart=qstart
}
' "$SORTED" > "$GROUPED"

# ------------------------------
# Calculate query length
# ------------------------------
FASTA_FILE="../AGO2.fasta"
if [[ -f "$FASTA_FILE" ]]; then
    QUERY_LENGTH=$(grep -v "^>" "$FASTA_FILE" | tr -d '\n' | wc -c)
else
    echo "Warning: $FASTA_FILE not found"
    QUERY_LENGTH=0
fi
echo "AGO2 query length: $QUERY_LENGTH"

# ------------------------------
# Weighted percent identity
# ------------------------------
IDENTITY_FILE="${BASENAME}_combined.tsv.sorted.grouped.identity"
echo "Calculating weighted percent identity..."
awk -F'\t' '
NR==1 { print $0 "\thit_percent_identity"; next }
{
    hit_id=$3; len=$9+0; percent_id=$12+0
    sum_len[hit_id]+=len
    sum_weighted[hit_id]+=len*percent_id
    rows[NR]=$0
    hit_ids[NR]=hit_id
    max_line=NR
}
END {
    for(h in sum_len) weighted[h]=sum_weighted[h]/sum_len[h]
    for(i=2;i<=max_line;i++) printf "%s\t%.1f\n", rows[i], weighted[hit_ids[i]]
}
' "$GROUPED" > "$IDENTITY_FILE"

# ------------------------------
# Coverage calculation
# ------------------------------
COVERAGE_FILE="${IDENTITY_FILE}.coverage"
echo "Calculating coverage..."
awk -v qlen="$QUERY_LENGTH" -F'\t' -v OFS='\t' '
function merge_intervals(arr, n, merged, m, i, cur_start, cur_end, s, e) {
    m=0; cur_start=-1; cur_end=-1
    for(i=1;i<=n;i++){
        split(arr[i], parts, ","); s=parts[1]; e=parts[2]
        if(cur_start<0){ cur_start=s; cur_end=e }
        else if(s<=cur_end){ if(e>cur_end) cur_end=e }
        else { merged[++m]=cur_start "," cur_end; cur_start=s; cur_end=e }
    }
    merged[++m]=cur_start "," cur_end
    return m
}
NR==1 { print $0, "hit_coverage_percent"; next }
{
    hit_id=$3; qs=$5+0; qe=$6+0
    if(qs>qe){ tmp=qs; qs=qe; qe=tmp }
    count[hit_id]++; intervals[hit_id "," count[hit_id]]=qs "," qe
    lines[NR]=$0; line_hit[NR]=hit_id; max_line=NR
}
END {
    for(h in count){
        n=count[h]
        for(i=1;i<=n;i++){ split(intervals[h","i],p,","); start[i]=p[1]; end[i]=p[2] }
        # sort intervals by start
        for(i=1;i<=n;i++){ for(j=i+1;j<=n;j++){ if(start[i]>start[j]){ tmp=start[i]; start[i]=start[j]; start[j]=tmp; tmp=end[i]; end[i]=end[j]; end[j]=tmp } } }
        for(i=1;i<=n;i++) sorted[i]=start[i] "," end[i]
        m=merge_intervals(sorted,n,merged)
        cov=0
        for(i=1;i<=m;i++){ split(merged[i],p,","); cov+=(p[2]-p[1]+1) }
        coverage[h]=(cov/qlen)*100
        delete start; delete end; delete sorted; delete merged
    }
    for(i=2;i<=max_line;i++){
        h=line_hit[i]
        printf "%s\t%.1f\n", lines[i], (coverage[h]?coverage[h]:0)
    }
}
' "$IDENTITY_FILE" > "$COVERAGE_FILE"

echo "Done! Output files:"
echo "  - $SORTED"
echo "  - $GROUPED"
echo "  - $IDENTITY_FILE"
echo "  - $COVERAGE_FILE"
