from pathlib import Path

configfile: "config/config.yaml"

GENE = config["gene"]
OUTDIR = Path(config.get("output_dir", "results")) / GENE
CLASSIFIER = config.get("classifier_script", "Scripts_And_Guide/classify_hits_v2.py")


def all_targets():
    targets = [
        OUTDIR / "01_tblastn" / f"{GENE}_combined.tsv.sorted.grouped.identity.coverage",
        OUTDIR / "02_aa_sequence" / f"{GENE}_newhits_filtered.txt",
        OUTDIR / "02_aa_sequence" / f"{GENE}_newhits_filtered.fasta",
        OUTDIR / "03_reciprocal_tblastn" / f"{GENE}_newhits_filtered_confirmed.fasta",
        OUTDIR / "04_domain_analysis" / f"{GENE}_classification.tsv",
        OUTDIR / "05_check_contig_ends" / f"{GENE}_partial_copies_annotated_summary.tsv",
        OUTDIR / "06_predict_unannotated" / f"{GENE}_predicted_proteins_clean.fa",
        OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed.fasta",
        OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed_partial_summary.tsv",
    ]
    if config.get("run_missing_domains", False):
        targets.append(OUTDIR / "07_find_missing_domains" / f"{GENE}_record_keeping_missing_domains.txt")
    return targets

rule all:
    input:
        all_targets()

rule tblastn_and_percentages:
    input:
        query=lambda wc: config["query_fasta"]
    output:
        combined=OUTDIR / "01_tblastn" / f"{GENE}_combined.tsv",
        sorted=OUTDIR / "01_tblastn" / f"{GENE}_combined.tsv.sorted",
        grouped=OUTDIR / "01_tblastn" / f"{GENE}_combined.tsv.sorted.grouped",
        identity=OUTDIR / "01_tblastn" / f"{GENE}_combined.tsv.sorted.grouped.identity",
        coverage=OUTDIR / "01_tblastn" / f"{GENE}_combined.tsv.sorted.grouped.identity.coverage"
    shell:
        "python3 workflow/scripts/tblastn_and_percentages.py "
        "--query {input.query} "
        "--gene {GENE} "
        "--genome-dir {config[references][genome_dir]} "
        "--genome-glob '{config[references][genome_glob]}' "
        "--combined {output.combined} "
        "--sorted {output.sorted} "
        "--grouped {output.grouped} "
        "--identity {output.identity} "
        "--coverage {output.coverage}"

rule filter_hits:
    input:
        OUTDIR / "01_tblastn" / f"{GENE}_combined.tsv.sorted.grouped.identity.coverage"
    output:
        OUTDIR / "02_aa_sequence" / f"{GENE}_newhits_filtered.txt"
    shell:
        "python3 workflow/scripts/filter_hits.py "
        "--input {input} "
        "--output {output} "
        "--identity-threshold {config[thresholds][identity]} "
        "--coverage-threshold {config[thresholds][coverage]}"

rule extract_annotated_proteins:
    input:
        hits=OUTDIR / "02_aa_sequence" / f"{GENE}_newhits_filtered.txt"
    output:
        fasta=OUTDIR / "02_aa_sequence" / f"{GENE}_newhits_filtered.fasta",
        nocds=OUTDIR / "02_aa_sequence" / f"{GENE}_newhits_filtered_noCDS.bed"
    shell:
        "python3 workflow/scripts/extract_annotated_proteins.py "
        "--input {input.hits} "
        "--output-fasta {output.fasta} "
        "--output-bed {output.nocds} "
        "--genome-dir {config[references][genome_dir]} "
        "--annotation-suffix {config[references][annotation_suffix]} "
        "--genome-suffix {config[references][genome_suffix]}"

rule reciprocal_annotated:
    input:
        fasta=OUTDIR / "02_aa_sequence" / f"{GENE}_newhits_filtered.fasta"
    output:
        confirmed=OUTDIR / "03_reciprocal_tblastn" / f"{GENE}_newhits_filtered_confirmed.fasta",
        eliminated=OUTDIR / "03_reciprocal_tblastn" / f"{GENE}_newhits_filtered_eliminated.fasta"
    shell:
        "python3 workflow/scripts/reciprocal_blast.py "
        "--input {input.fasta} "
        "--confirmed {output.confirmed} "
        "--eliminated {output.eliminated} "
        "--reference-fasta {config[reciprocal][fasta]} "
        "--blast-db {config[reciprocal][db]} "
        "--gff {config[reciprocal][gff]} "
        "--gene-match {config[reciprocal][gene_match]}"

rule hmmscan_annotated:
    input:
        fasta=OUTDIR / "03_reciprocal_tblastn" / f"{GENE}_newhits_filtered_confirmed.fasta",
        species_list=lambda wc: config["species_list"]
    output:
        domtbl=OUTDIR / "04_domain_analysis" / f"{GENE}_domtblout.tsv",
        classification=OUTDIR / "04_domain_analysis" / f"{GENE}_classification.tsv",
        species_summary=OUTDIR / "04_domain_analysis" / f"{GENE}_species_summary.tsv",
        full=OUTDIR / "04_domain_analysis" / f"{GENE}_full_copies.fasta",
        partial=OUTDIR / "04_domain_analysis" / f"{GENE}_partial_copies.fasta",
        ignore=OUTDIR / "04_domain_analysis" / f"{GENE}_ignore_copies.fasta"
    shell:
        "python3 workflow/scripts/run_hmmscan_and_split.py "
        "--gene {GENE} "
        "--input {input.fasta} "
        "--pfam {config[references][pfam_db]} "
        "--domtbl {output.domtbl} "
        "--classification {output.classification} "
        "--species-summary {output.species_summary} "
        "--species-list {input.species_list} "
        "--classifier {CLASSIFIER} "
        "--full {output.full} "
        "--partial {output.partial} "
        "--ignore {output.ignore}"

rule longest_introns:
    input:
        OUTDIR / "04_domain_analysis" / f"{GENE}_full_copies.fasta"
    output:
        OUTDIR / "05_check_contig_ends" / f"{GENE}_longest_introns.txt"
    shell:
        "python3 workflow/scripts/compute_intron_threshold.py "
        "--input {input} "
        "--output {output} "
        "--genome-dir {config[references][genome_dir]} "
        "--annotation-suffix {config[references][annotation_suffix]}"

rule check_contig_ends_annotated:
    input:
        partial=OUTDIR / "04_domain_analysis" / f"{GENE}_partial_copies.fasta",
        threshold=OUTDIR / "05_check_contig_ends" / f"{GENE}_longest_introns.txt"
    output:
        likely=OUTDIR / "05_check_contig_ends" / f"{GENE}_partial_copies_annotated_likely_full.fasta",
        remaining=OUTDIR / "05_check_contig_ends" / f"{GENE}_partial_copies_annotated_remaining.fasta",
        summary=OUTDIR / "05_check_contig_ends" / f"{GENE}_partial_copies_annotated_summary.tsv"
    shell:
        "python3 workflow/scripts/check_contig_ends.py "
        "--mode annotated "
        "--input {input.partial} "
        "--threshold-file {input.threshold} "
        "--likely-full {output.likely} "
        "--remaining {output.remaining} "
        "--summary {output.summary} "
        "--genome-dir {config[references][genome_dir]} "
        "--annotation-suffix {config[references][annotation_suffix]} "
        "--genome-suffix {config[references][genome_suffix]}"

rule augustus_predict:
    input:
        bed=OUTDIR / "02_aa_sequence" / f"{GENE}_newhits_filtered_noCDS.bed"
    output:
        proteins=OUTDIR / "06_predict_unannotated" / f"{GENE}_predicted_proteins_clean.fa"
    shell:
        "python3 workflow/scripts/run_augustus_extract.py "
        "--input-bed {input.bed} "
        "--output-fasta {output.proteins} "
        "--genome-dir {config[references][genome_dir]} "
        "--genome-suffix {config[references][genome_suffix]} "
        "--augustus-species {config[augustus][species]} "
        "--extend-bp {config[augustus][extend_bp]}"

rule reciprocal_predicted:
    input:
        fasta=OUTDIR / "06_predict_unannotated" / f"{GENE}_predicted_proteins_clean.fa"
    output:
        confirmed=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed.fasta",
        eliminated=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_eliminated.fasta"
    shell:
        "python3 workflow/scripts/reciprocal_blast.py "
        "--input {input.fasta} "
        "--confirmed {output.confirmed} "
        "--eliminated {output.eliminated} "
        "--reference-fasta {config[reciprocal][fasta]} "
        "--blast-db {config[reciprocal][db]} "
        "--gff {config[reciprocal][gff]} "
        "--gene-match {config[reciprocal][gene_match]}"

rule hmmscan_predicted:
    input:
        fasta=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed.fasta",
        species_list=lambda wc: config["species_list"]
    output:
        domtbl=OUTDIR / "06_predict_unannotated" / f"{GENE}_predicted_domtblout.tsv",
        classification=OUTDIR / "06_predict_unannotated" / f"{GENE}_predicted_classification.tsv",
        species_summary=OUTDIR / "06_predict_unannotated" / f"{GENE}_predicted_species_summary.tsv",
        full=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed_full_copies.fasta",
        partial=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed_partial_copies.fasta",
        ignore=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed_ignore_copies.fasta"
    shell:
        "python3 workflow/scripts/run_hmmscan_and_split.py "
        "--gene {GENE} "
        "--input {input.fasta} "
        "--pfam {config[references][pfam_db]} "
        "--domtbl {output.domtbl} "
        "--classification {output.classification} "
        "--species-summary {output.species_summary} "
        "--species-list {input.species_list} "
        "--classifier {CLASSIFIER} "
        "--full {output.full} "
        "--partial {output.partial} "
        "--ignore {output.ignore}"

rule check_contig_ends_predicted:
    input:
        partial=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed_partial_copies.fasta",
        threshold=OUTDIR / "05_check_contig_ends" / f"{GENE}_longest_introns.txt"
    output:
        likely=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed_partial_likelyfull.fasta",
        remaining=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed_partial_remaining.fasta",
        summary=OUTDIR / "06_predict_unannotated" / f"{GENE}_NOCDS_AUGUSTUS_confirmed_partial_summary.tsv"
    shell:
        "python3 workflow/scripts/check_contig_ends.py "
        "--mode predicted "
        "--input {input.partial} "
        "--threshold-file {input.threshold} "
        "--likely-full {output.likely} "
        "--remaining {output.remaining} "
        "--summary {output.summary} "
        "--genome-dir {config[references][genome_dir]} "
        "--annotation-suffix {config[references][annotation_suffix]} "
        "--genome-suffix {config[references][genome_suffix]}"

rule missing_domains:
    input:
        table=lambda wc: config["missing_domains"]["input_table"]
    output:
        OUTDIR / "07_find_missing_domains" / f"{GENE}_record_keeping_missing_domains.txt"
    params:
        strand_column=lambda wc: config["missing_domains"].get("strand_column", "STRAND")
    shell:
        "python3 workflow/scripts/find_missing_domains_strand_aware.py "
        "--gene {GENE} "
        "--input-table {input.table} "
        "--output-table {output} "
        "--pfam {config[references][pfam_db]} "
        "--genome-dir {config[references][genome_dir]} "
        "--genome-suffix {config[references][genome_suffix]} "
        "--strand-column {params.strand_column}"
