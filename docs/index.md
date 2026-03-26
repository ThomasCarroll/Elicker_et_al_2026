---
title: Pipeline Overview
nav_order: 1
---

# Fly Pipeline

This workflow packages the fly gene-copy pipeline as a Dockerized Snakemake project. The aim is to make the original process reproducible, configurable, and easier to rerun without editing shell scripts by hand for each gene or dataset.

## What the pipeline does

The workflow starts from a protein query FASTA and searches a set of genomes for homologous regions. It then filters candidate hits, extracts overlapping annotated protein-coding sequences where possible, validates candidates by reciprocal search against a reference genome, classifies them by expected protein domain architecture, checks whether partial copies sit near contig ends, and optionally predicts unannotated proteins in no-CDS regions with AUGUSTUS.

## Workflow stages

### 1. Search genomes with `tblastn`

- Runs `tblastn` against each configured genome
- Collects hits into a combined table
- Computes weighted percent identity and merged query coverage per grouped hit

### 2. Filter and extract annotated proteins

- Applies identity and coverage thresholds
- Finds overlapping annotated transcripts and CDS features
- Translates extracted CDS sequences into protein FASTA
- Records regions with no overlapping CDS as candidates for prediction

### 3. Reciprocal validation of annotated hits

- Searches extracted proteins back against the reciprocal reference genome
- Splits sequences into confirmed and eliminated sets based on overlap with the target gene in the reciprocal annotation

### 4. Domain analysis

- Runs `hmmscan` against Pfam
- Classifies sequences as `FULL`, `PARTIAL`, or `IGNORE`
- Splits FASTA outputs by class

### 5. Contig-end checks for annotated partials

- Computes an intron-length threshold from the confirmed full copies
- Marks partial copies near contig ends as likely full
- Produces a summary table for the annotated branch

### 6. Predict proteins in no-CDS regions

- Extracts genomic sequence around regions with no overlapping CDS
- Runs AUGUSTUS on those regions
- Validates predicted proteins with reciprocal search and domain analysis
- Applies the same contig-end logic to predicted partials

### 7. Optional missing-domain rescue

- Re-examines partial hits by scanning an extended genomic region around them
- Uses a strand-aware three-frame translation rather than all six frames
- Reports whether missing domains support reclassification

## Inputs

The workflow expects:

- a query FASTA for the target gene
- genome FASTA files
- matching GFF annotation files
- a reciprocal reference genome plus annotation
- a Pfam database and its HMMER index files

## Outputs

Outputs are written under `results/<GENE>/` by stage. The main output families are:

- combined search tables
- filtered hit tables
- annotated protein FASTA files
- confirmed and eliminated reciprocal hit FASTA files
- domain classification tables
- contig-end summary tables
- AUGUSTUS-predicted protein FASTA files

## Project layout

- `Snakefile`: workflow entry point
- `config/`: configuration files
- `workflow/scripts/`: parameterized helper scripts
- `tests/`: smoke and full workflow tests
- `Scripts_And_Guide/`: original scripts and guide kept for reference
