#!/usr/bin/env bash
set -euo pipefail

IMAGE_NAME="${1:-kimberly}"
WORKDIR="/work"
CONFIG_FILE="${CONFIG_FILE:-config/config.yaml}"

read_config_value() {
  python3 - "$1" "$2" <<'PY'
import yaml
import sys
with open(sys.argv[1]) as fh:
    cfg = yaml.safe_load(fh)
print(cfg[sys.argv[2]])
PY
}

GENE="${GENE:-$(read_config_value "${CONFIG_FILE}" gene)}"

run_in_container() {
  docker run --rm -t \
    -v "$PWD:${WORKDIR}" \
    -w "${WORKDIR}" \
    "${IMAGE_NAME}" "$@"
}

check_file() {
  [[ -f "$1" ]] || {
    echo "Expected output not found: $1" >&2
    exit 1
  }
}

echo "Dry run"
run_in_container snakemake --configfile "${CONFIG_FILE}" -n -p

echo "Stage 1"
run_in_container snakemake --configfile "${CONFIG_FILE}" -j 1 -p "results/${GENE}/01_tblastn/${GENE}_combined.tsv.sorted.grouped.identity.coverage"

echo "Stage 2"
run_in_container snakemake --configfile "${CONFIG_FILE}" -j 1 -p "results/${GENE}/02_aa_sequence/${GENE}_newhits_filtered.fasta"

echo "Stage 4"
run_in_container snakemake --configfile "${CONFIG_FILE}" -j 1 -p "results/${GENE}/04_domain_analysis/${GENE}_classification.tsv"

check_file "results/${GENE}/01_tblastn/${GENE}_combined.tsv.sorted.grouped.identity.coverage"
check_file "results/${GENE}/02_aa_sequence/${GENE}_newhits_filtered.fasta"
check_file "results/${GENE}/04_domain_analysis/${GENE}_classification.tsv"

echo "Smoke test passed"
