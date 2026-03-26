# Fly Pipeline

This repository contains a Dockerized Snakemake workflow for running the fly gene-copy pipeline from a single configuration file instead of editing the original shell scripts by hand for each run.

Original pipeline and analysis framework by Kimberly Elicker.

## Repository layout

- `Snakefile`: workflow entry point
- `config/`: run configuration and species list
- `docs/`: GitHub Pages content for overview, tutorial, and install/test docs
- `workflow/scripts/`: parameterized helper scripts used by Snakemake
- `tests/`: smoke and full workflow test scripts
- `Dockerfile` and `environment.yml`: container build and software environment
- `Scripts_And_Guide/`: original scripts and guide retained for provenance and comparison

## Tracked code vs local data

The repository is intended to track workflow code, configuration, and tests. Large reference inputs and generated outputs are intentionally excluded from Git:

- genome FASTA files
- GFF annotation files
- Pfam database files
- BLAST index files
- `results/`
- `.snakemake/`

That keeps the repository suitable for GitHub while still allowing the workflow to run locally against external reference data.

## Quick start

Build the image:

```bash
docker build -t fly-pipeline .
```

Run the workflow:

```bash
docker run --rm -it -v "$PWD:/work" -w /work fly-pipeline snakemake -j 4
```

The workflow reads its settings from `config/config.yaml`.

## Tests

Smoke test:

```bash
bash tests/smoke_test.sh fly-pipeline
```

Full test:

```bash
bash tests/full_test.sh fly-pipeline
```

## Notes

- The current example config uses `Drosophila_affinis_*` as the reciprocal reference because those files are present in this working directory.
- If you want to mirror the original guide exactly, update `config/config.yaml` to point at your intended reference species.
- The Docker image is built from the workflow code only. Reference data are expected to be present in the mounted working directory at runtime.
