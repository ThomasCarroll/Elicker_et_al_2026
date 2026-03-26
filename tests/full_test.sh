#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME="${1:-kimberly}"
WORKDIR="/work"
CONFIG_FILE="${CONFIG_FILE:-config/config.yaml}"

read_config() {
  python3 - "$1" "$2" <<'PY'
import yaml
import sys
with open(sys.argv[1]) as fh:
    cfg = yaml.safe_load(fh)
value = cfg
for key in sys.argv[2].split("."):
    value = value[key]
print(value)
PY
}

require_file() {
  if [[ ! -f "$1" ]]; then
    echo "Missing required file: $1" >&2
    exit 1
  fi
}

run_in_container() {
  docker run --rm -t \
    -v "$PWD:${WORKDIR}" \
    -w "${WORKDIR}" \
    "${IMAGE_NAME}" "$@"
}

require_file "${CONFIG_FILE}"
require_file "Snakefile"
require_file "Pfam-A.hmm"

GENE="$(read_config "${CONFIG_FILE}" gene)"
RUN_MISSING="$(read_config "${CONFIG_FILE}" run_missing_domains)"

echo "Dry run"
run_in_container snakemake --configfile "${CONFIG_FILE}" -n -p

echo "Full workflow"
run_in_container snakemake --configfile "${CONFIG_FILE}" -j "${JOBS:-2}" -p

require_file "results/${GENE}/01_tblastn/${GENE}_combined.tsv.sorted.grouped.identity.coverage"
require_file "results/${GENE}/02_aa_sequence/${GENE}_newhits_filtered.txt"
require_file "results/${GENE}/02_aa_sequence/${GENE}_newhits_filtered.fasta"
require_file "results/${GENE}/03_reciprocal_tblastn/${GENE}_newhits_filtered_confirmed.fasta"
require_file "results/${GENE}/04_domain_analysis/${GENE}_classification.tsv"
require_file "results/${GENE}/05_check_contig_ends/${GENE}_partial_copies_annotated_summary.tsv"
require_file "results/${GENE}/06_predict_unannotated/${GENE}_predicted_proteins_clean.fa"
require_file "results/${GENE}/06_predict_unannotated/${GENE}_NOCDS_AUGUSTUS_confirmed.fasta"
require_file "results/${GENE}/06_predict_unannotated/${GENE}_NOCDS_AUGUSTUS_confirmed_partial_summary.tsv"

if [[ "${RUN_MISSING}" == "True" || "${RUN_MISSING}" == "true" ]]; then
  require_file "results/${GENE}/07_find_missing_domains/${GENE}_record_keeping_missing_domains.txt"
fi

echo "Full test passed"
