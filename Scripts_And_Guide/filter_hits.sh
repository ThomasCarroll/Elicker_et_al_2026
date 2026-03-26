#!/bin/bash

gene_name="PIWI"

# Input file
input="${gene_name}_newhits.txt"

# Output file
output="${gene_name}_newhits_filtered.txt"

# Thresholds
identity_thresh=55       # minimum percent identity
coverage_thresh=20       # minimum percent coverage

awk -v gene="$gene_name" -v input="$input" -v id="$identity_thresh" -v cov="$coverage_thresh" '
BEGIN { 
    FS=OFS="\t"; 
    print "Starting filtering for gene:", gene > "/dev/stderr"
    print "Input file:", input > "/dev/stderr"
    print "Identity threshold:", id > "/dev/stderr"
    print "Coverage threshold:", cov > "/dev/stderr"
}
NR==1 { print; next }  # print header

{
    # Keep track of all unique hit_ids
    all_hits[$3]=1

    # Check thresholds
    if ($16 > id && $17 > cov) {
        passed_hits[$3]=1
        print
    }
}

END {
    print "Total unique hit_ids:", length(all_hits) > "/dev/stderr"
    print "Unique hit_ids passing thresholds:", length(passed_hits) > "/dev/stderr"
}' "$input" > "$output"

echo "Filtering complete. Filtered hits saved to $output"
