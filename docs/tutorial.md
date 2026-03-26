---
title: Tutorial
nav_order: 2
---

# Tutorial

This page walks through the normal way to run the workflow on a local checkout.

## 1. Prepare the working directory

The workflow code is meant to live in Git, but the large reference files should stay local. Make sure the following are present in your working directory:

- genome FASTA files
- matching GFF files
- the Pfam HMM database
- the reciprocal reference genome and annotation

The example configuration in this repository uses `Drosophila_affinis_*` for the reciprocal step because those files are present in the test dataset.

## 2. Check the configuration

Edit `config/config.yaml` and review:

- `gene`
- `query_fasta`
- `references.*`
- `reciprocal.*`
- `thresholds.*`
- `augustus.*`

In particular, update the reciprocal reference if you want to mirror the original guide with a different species.

## 3. Build the Docker image

```bash
docker build -t fly-pipeline .
```

## 4. Run a dry run

Use Snakemake's dry-run mode first to verify the DAG and file paths:

```bash
docker run --rm -it -v "$PWD:/work" -w /work fly-pipeline snakemake -n -p
```

## 5. Run the workflow

```bash
docker run --rm -it -v "$PWD:/work" -w /work fly-pipeline snakemake -j 4
```

You can lower or raise `-j` depending on the machine.

## 6. Inspect outputs

The workflow writes outputs under:

```text
results/<GENE>/
```

Useful places to inspect first:

- `01_tblastn/`
- `02_aa_sequence/`
- `03_reciprocal_tblastn/`
- `04_domain_analysis/`
- `05_check_contig_ends/`
- `06_predict_unannotated/`

## 7. Re-run after a config change

If you change parameters or inputs, rerun Snakemake from the same command. It will only rebuild what is out of date.

To force a clean rerun:

```bash
rm -rf results
docker run --rm -it -v "$PWD:/work" -w /work fly-pipeline snakemake -j 4
```

## Notes on empty intermediate files

Some branches can legitimately produce no candidates. For example, reciprocal validation may yield no confirmed hits for a dataset. The workflow is set up to handle empty downstream inputs without crashing.
