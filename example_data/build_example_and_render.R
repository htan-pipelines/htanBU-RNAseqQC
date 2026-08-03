#!/usr/bin/env Rscript
# ==============================================================================
# build_example_and_render.R
#
# Builds a small SYNTHETIC dataset that mimics the four inputs qc_report.Rmd
# expects (aggregated SummarizedExperiment, external sample metadata, somalier
# pairwise relatedness, arcasHLA genotypes), writes them to example_data/ as
# real files (.rds/.tsv, read back the same way render_qc_report.R does), and
# renders qc_report.Rmd against them end to end. Used as a smoke test -- not
# real biological data.
#
# Usage:
#   Rscript example_data/build_example_and_render.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(rmarkdown)
})

set.seed(42)

repo_dir <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))), ".."))
if (!file.exists(file.path(repo_dir, "qc_report.Rmd"))) repo_dir <- getwd()  # fallback for interactive use
out_dir <- file.path(repo_dir, "example_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Sample design
# ------------------------------------------------------------------------------
# qc_report.Rmd's se.subset filters by tissue_type only, so all "Lung" samples
# (across both sites below) end up in the tissue we'll actually render;
# collection_site still varies within/across tissues purely to give the general
# stats plots and PCA shape aesthetics something to show.
n_lung_siteA <- 10
n_lung_siteB <- 2
n_nasal_siteA <- 2
n <- n_lung_siteA + n_lung_siteB + n_nasal_siteA

patient_id <- sprintf("SynCohort_P%02d", seq_len(n))
sample_id  <- paste0(patient_id, "_S1")

tissue_type <- c(rep("Lung", n_lung_siteA + n_lung_siteB), rep("Nasal", n_nasal_siteA))
collection_site <- c(rep("SiteA", n_lung_siteA), rep("SiteB", n_lung_siteB), rep("SiteA", n_nasal_siteA))
cohort <- rep(c("CohortX", "CohortY"), length.out = n)
batch_id <- rep(c("Batch1", "Batch2"), length.out = n)

# Balanced sex assignment (used to build the sex-marker expression signal below)
sex <- rep(c("Female", "Male"), length.out = n)

rin <- round(runif(n, 4, 9.5), 1)
dv200 <- round(runif(n, 45, 92), 1)
tin_median <- round(runif(n, 40, 92), 1)

# ------------------------------------------------------------------------------
# 2. colData QC columns (fastqc_/STAR_/rnaseqc_/samtools_/TIN_-prefixed, as
#    produced by RNA_seq_pipeline.wdl + combine_se.wdl) -- everything
#    qc_report.Rmd reads off colData(se) by exact name.
# ------------------------------------------------------------------------------
total_reads <- round(runif(n, 2e7, 6e7))
mapped_reads <- round(total_reads * runif(n, 0.85, 0.98))
mapped_unique_reads <- round(mapped_reads * runif(n, 0.7, 0.95))
high_quality_reads <- round(total_reads * runif(n, 0.8, 0.97))
rrna_reads <- round(total_reads * runif(n, 0.001, 0.02))
star_uniquely_mapped <- round(total_reads * runif(n, 0.75, 0.95))
star_multimapped_multiple <- round(total_reads * runif(n, 0.01, 0.05))

qc_cols <- data.frame(
  fastqc_per_base_seq_quality        = sample(c("PASS", "WARN"), n, replace = TRUE),
  Somalier_n_high_relatedness_pairs  = 0L,
  Genotypes_match_mean               = round(runif(n, 0, 0.3), 2),
  samtools_reads_properly_paired     = round(mapped_reads * runif(n, 0.9, 0.99)),
  TIN_median                         = tin_median,
  rnaseqc_Sample                     = sample_id,
  rnaseqc_Total.Reads                = total_reads,
  rnaseqc_Mapped.Reads                = mapped_reads,
  rnaseqc_Mapped.Unique.Reads         = mapped_unique_reads,
  rnaseqc_Unique.Rate.of.Mapped       = round(mapped_unique_reads / mapped_reads, 3),
  rnaseqc_High.Quality.Reads          = high_quality_reads,
  rnaseqc_Low.Quality.Reads           = total_reads - high_quality_reads,
  rnaseqc_rRNA.Reads                  = rrna_reads,
  rnaseqc_rRNA.Rate                   = round(rrna_reads / total_reads, 4),
  rnaseqc_Median.3..bias              = round(runif(n, 0.2, 0.7), 3),
  rnaseqc_Genes.Detected              = round(runif(n, 8000, 18000)),
  rnaseqc_Median.Exon.CV              = round(runif(n, 0.3, 1.3), 3),
  rnaseqc_Duplicate.Rate.of.Mapped    = round(runif(n, 0.1, 0.5), 3),
  rnaseqc_Median.of.Avg.Transcript.Coverage = round(runif(n, 5, 40), 2),
  STAR_total_reads                    = total_reads,
  STAR_uniquely_mapped                = star_uniquely_mapped,
  STAR_uniquely_mapped_percent        = round(star_uniquely_mapped / total_reads * 100, 2),
  STAR_multimapped_multiple           = star_multimapped_multiple,
  STAR_multimapped_multiple_percent   = round(star_multimapped_multiple / total_reads * 100, 2),
  STAR_unmapped_tooshort_percent      = round(runif(n, 0.5, 5), 2),
  STAR_avg_mapped_read_length         = round(runif(n, 90, 150), 1),
  STAR_num_splices                    = round(runif(n, 1e6, 5e6)),
  row.names = sample_id,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 3. Genes: 7 sex-marker control genes (XIST, and the 6 Y-linked genes the
#    report looks for) + ~35 background genes -- "a few dozen" total.
# ------------------------------------------------------------------------------
sex_genes <- c("XIST", "TXLNGY", "DDX3Y", "KDM5D", "RPS4Y1", "USP9Y", "UTY")
bg_genes <- paste0("GENE", seq_len(35))
gene_symbols <- c(sex_genes, bg_genes)
n_genes <- length(gene_symbols)

row_data <- DataFrame(
  gene_id            = sprintf("ENSG%011d", seq_len(n_genes)),
  external_gene_name = gene_symbols,
  hgnc_symbol        = gene_symbols,
  entrezgene_id      = seq_len(n_genes),
  transcript_count   = sample(1:6, n_genes, replace = TRUE),
  gene_biotype       = "protein_coding",
  chromosome_name    = ifelse(gene_symbols == "XIST", "chrX",
                        ifelse(gene_symbols %in% sex_genes, "chrY",
                               sample(paste0("chr", 1:22), n_genes, replace = TRUE))),
  gene_length        = sample(500:8000, n_genes, replace = TRUE),
  row.names           = gene_symbols
)

# ------------------------------------------------------------------------------
# 4. Expression counts: background genes get sample-varying negative-binomial
#    counts; sex genes get a strong XIST-vs-Y-linked signal split by `sex` so
#    the sex-check section (kmeans on PC1 of sex-marker expression) has a real
#    bimodal signal to find, instead of being pure noise.
# ------------------------------------------------------------------------------
lib_factor <- runif(n, 0.7, 1.4)  # per-sample library-size-ish scaling

bg_counts <- t(vapply(bg_genes, function(g) {
  mu <- runif(1, 50, 4000) * lib_factor
  rnbinom(n, mu = mu, size = 8)
}, numeric(n)))
rownames(bg_counts) <- bg_genes

sex_counts <- matrix(0, nrow = length(sex_genes), ncol = n, dimnames = list(sex_genes, sample_id))
is_female <- sex == "Female"
sex_counts["XIST", is_female]  <- rnbinom(sum(is_female),  mu = 900 * lib_factor[is_female],  size = 10)
sex_counts["XIST", !is_female] <- rnbinom(sum(!is_female), mu = 8   * lib_factor[!is_female],  size = 10)
for (g in setdiff(sex_genes, "XIST")) {
  sex_counts[g, is_female]  <- rnbinom(sum(is_female),  mu = 5   * lib_factor[is_female],  size = 10)
  sex_counts[g, !is_female] <- rnbinom(sum(!is_female), mu = 400 * lib_factor[!is_female], size = 10)
}

expected_count <- rbind(sex_counts, bg_counts)[gene_symbols, ]
storage.mode(expected_count) <- "double"
colnames(expected_count) <- sample_id

tpm <- sweep(expected_count, 2, colSums(expected_count), FUN = function(x, s) x / s * 1e6)
fpkm <- tpm  # not exact FPKM math, but qc_report.Rmd never reads this assay -- placeholder is fine

se <- SummarizedExperiment(
  assays = list(expected_count = expected_count, TPM = tpm, FPKM = fpkm),
  rowData = row_data,
  colData = DataFrame(qc_cols)
)

# ------------------------------------------------------------------------------
# 5. External per-sample metadata (rnaAnnot) -- column names match the
#    defaults in render_qc_report.R's column_map, so no remapping needed here.
# ------------------------------------------------------------------------------
rnaAnnot <- data.frame(
  sample_id = sample_id,
  patient_id = patient_id,
  tissue_type = tissue_type,
  collection_site = collection_site,
  cohort = cohort,
  batch_id = batch_id,
  rin = rin,
  dv200 = dv200,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 6. Somalier pairwise relatedness (somalier.pairs.tsv). Real somalier output
#    has one row per unordered sample pair, with a leading "#sample_a" header
#    that read.delim() turns into "X.sample_a" -- write it out the same way
#    and read it back, instead of constructing the post-read.delim column name
#    directly, so this exercises the actual file-reading path.
# ------------------------------------------------------------------------------
pairs <- combn(sample_id, 2)
n_pairs <- ncol(pairs)
hets_a <- round(runif(n_pairs, 2000, 5000))
hom_alts_a <- round(hets_a * runif(n_pairs, 0.3, 0.7))
hets_b <- round(runif(n_pairs, 2000, 5000))
hom_alts_b <- round(hets_b * runif(n_pairs, 0.3, 0.7))

somalier_pairs_out <- data.frame(
  "#sample_a" = pairs[1, ],
  sample_b = pairs[2, ],
  relatedness = round(rbeta(n_pairs, 1.5, 8), 3),   # concentrated well below the 0.6 swap cutoff
  hets_a = hets_a,
  hom_alts_a = hom_alts_a,
  hets_b = hets_b,
  hom_alts_b = hom_alts_b,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 7. arcasHLA genotype calls (genotypes.tsv): "subject" + a 7-gene x 2-allele
#    panel (14 columns), matching the width the report's HLA match score
#    (shared_features / ncol(x)) is designed to adapt to.
# ------------------------------------------------------------------------------
hla_loci <- c("A", "B", "C", "DPB1", "DQA1", "DQB1", "DRB1")
allele_pool <- function(locus) paste0(locus, "*", sprintf("%02d:%02d", sample(1:20, 6), sample(1:99, 6)))

genotypes_out <- data.frame(subject = sample_id, stringsAsFactors = FALSE)
for (locus in hla_loci) {
  pool <- allele_pool(locus)
  genotypes_out[[paste0(locus, "_1")]] <- sample(pool, n, replace = TRUE)
  genotypes_out[[paste0(locus, "_2")]] <- sample(pool, n, replace = TRUE)
}

# ------------------------------------------------------------------------------
# 8. Write all four inputs to disk, then read them back exactly the way
#    render_qc_report.R does, so this is a genuine end-to-end file-based test.
# ------------------------------------------------------------------------------
se_path <- file.path(out_dir, "Aggregated_Gene_Expression.rds")
rnaAnnot_path <- file.path(out_dir, "metadata.tsv")
somalier_pairs_path <- file.path(out_dir, "somalier.pairs.tsv")
genotypes_path <- file.path(out_dir, "genotypes.tsv")

saveRDS(se, se_path)
write.table(rnaAnnot, rnaAnnot_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(somalier_pairs_out, somalier_pairs_path, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(genotypes_out, genotypes_path, sep = "\t", row.names = FALSE, quote = FALSE)

se_in <- readRDS(se_path)
rnaAnnot_in <- read.table(rnaAnnot_path, sep = "\t", header = TRUE)
somalier_pairs_in <- read.delim(somalier_pairs_path, header = TRUE, stringsAsFactors = FALSE)
genotypes_in <- read.delim(genotypes_path)

# ------------------------------------------------------------------------------
# 9. Render qc_report.Rmd against the synthetic inputs
# ------------------------------------------------------------------------------
output_dir <- file.path(out_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# qc_report.Rmd's last chunk does `rm(list = ls())` for a clean sessionInfo() --
# render() defaults to knitting in this script's own global env (envir =
# parent.frame()), so without an isolated envir that rm() would wipe this script's
# own variables out from under it. render_qc_report.R doesn't hit this because it
# calls render() from inside a function (a fresh frame each time); this top-level
# script needs to ask for isolation explicitly.
rmarkdown::render(
  file.path(repo_dir, "qc_report.Rmd"),
  envir = new.env(parent = globalenv()),
  params = list(
    se = se_in,
    rnaAnnot = rnaAnnot_in,
    somalier_pairs = somalier_pairs_in,
    genotypes = genotypes_in,
    column_map = list(
      sample_id = "sample_id", patient_id = "patient_id", tissue_type = "tissue_type",
      collection_site = "collection_site", cohort = "cohort", batch_id = "batch_id",
      rin = "rin", dv200 = "dv200"
    ),
    tissueType = "Synthetic Example QC",
    tissue = "Lung",
    output_se = file.path(output_dir, "annotated_Gene_Expression.rds"),
    output_se_per = file.path(output_dir, "Lung_Gene_Expression.rds"),
    qcFile = file.path(output_dir, "qc_flag_summary.tsv"),
    showSession = TRUE
  ),
  output_file = "example_QCReport.html",
  output_dir = output_dir,
  intermediates_dir = output_dir
)

message("Done. Report written to: ", file.path(output_dir, "example_QCReport.html"))

# ------------------------------------------------------------------------------
# 10. Also generate the standalone somalier relatedness network, same as
#     render_qc_report.R does after rendering the report.
# ------------------------------------------------------------------------------
source(file.path(repo_dir, "generate_somalier_network.R"))
generate_somalier_network(
  somalier_pairs = somalier_pairs_in,
  output_html = file.path(output_dir, "somalier_network_graph.html"),
  cutoff = 0.6
)

# ==============================================================================
# 11. A SECOND synthetic dataset, this time mimicking a cohort with NO
#     tissue_type/collection_site/cohort columns at all -- a real scenario
#     (a single-tissue study where the metadata sheet just doesn't track those
#     concepts) that exercises params$group_vars and the optional tissue_type
#     path end to end, on top of (not instead of) the dataset above.
# ==============================================================================
n2 <- 12
patient_id2 <- sprintf("Data2Cohort_P%02d", seq_len(n2))
sample_id2  <- paste0(patient_id2, "_S1")

# The study's actual variable of interest -- takes the place collection_site/cohort
# would otherwise play, via params$group_vars.
adequacy_group <- rep(c("Adequate", "Adequate", "Borderline", "Intermediate"), length.out = n2)
batch2 <- rep(c("Batch1", "Batch2"), length.out = n2)
sex2 <- rep(c("Female", "Male"), length.out = n2)

rin2 <- round(runif(n2, 4, 9.5), 1)
dv200_2 <- round(runif(n2, 45, 92), 1)
tin_median2 <- round(runif(n2, 40, 92), 1)

total_reads2 <- round(runif(n2, 2e7, 6e7))
mapped_reads2 <- round(total_reads2 * runif(n2, 0.85, 0.98))
mapped_unique_reads2 <- round(mapped_reads2 * runif(n2, 0.7, 0.95))
high_quality_reads2 <- round(total_reads2 * runif(n2, 0.8, 0.97))
rrna_reads2 <- round(total_reads2 * runif(n2, 0.001, 0.02))
star_uniquely_mapped2 <- round(total_reads2 * runif(n2, 0.75, 0.95))
star_multimapped_multiple2 <- round(total_reads2 * runif(n2, 0.01, 0.05))

qc_cols2 <- data.frame(
  fastqc_per_base_seq_quality        = sample(c("PASS", "WARN"), n2, replace = TRUE),
  Somalier_n_high_relatedness_pairs  = 0L,
  Genotypes_match_mean               = round(runif(n2, 0, 0.3), 2),
  samtools_reads_properly_paired     = round(mapped_reads2 * runif(n2, 0.9, 0.99)),
  TIN_median                         = tin_median2,
  rnaseqc_Sample                     = sample_id2,
  rnaseqc_Total.Reads                = total_reads2,
  rnaseqc_Mapped.Reads                = mapped_reads2,
  rnaseqc_Mapped.Unique.Reads         = mapped_unique_reads2,
  rnaseqc_Unique.Rate.of.Mapped       = round(mapped_unique_reads2 / mapped_reads2, 3),
  rnaseqc_High.Quality.Reads          = high_quality_reads2,
  rnaseqc_Low.Quality.Reads           = total_reads2 - high_quality_reads2,
  rnaseqc_rRNA.Reads                  = rrna_reads2,
  rnaseqc_rRNA.Rate                   = round(rrna_reads2 / total_reads2, 4),
  rnaseqc_Median.3..bias              = round(runif(n2, 0.2, 0.7), 3),
  rnaseqc_Genes.Detected              = round(runif(n2, 8000, 18000)),
  rnaseqc_Median.Exon.CV              = round(runif(n2, 0.3, 1.3), 3),
  rnaseqc_Duplicate.Rate.of.Mapped    = round(runif(n2, 0.1, 0.5), 3),
  rnaseqc_Median.of.Avg.Transcript.Coverage = round(runif(n2, 5, 40), 2),
  STAR_total_reads                    = total_reads2,
  STAR_uniquely_mapped                = star_uniquely_mapped2,
  STAR_uniquely_mapped_percent        = round(star_uniquely_mapped2 / total_reads2 * 100, 2),
  STAR_multimapped_multiple           = star_multimapped_multiple2,
  STAR_multimapped_multiple_percent   = round(star_multimapped_multiple2 / total_reads2 * 100, 2),
  STAR_unmapped_tooshort_percent      = round(runif(n2, 0.5, 5), 2),
  STAR_avg_mapped_read_length         = round(runif(n2, 90, 150), 1),
  STAR_num_splices                    = round(runif(n2, 1e6, 5e6)),
  row.names = sample_id2,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

lib_factor2 <- runif(n2, 0.7, 1.4)

bg_counts2 <- t(vapply(bg_genes, function(g) {
  mu <- runif(1, 50, 4000) * lib_factor2
  rnbinom(n2, mu = mu, size = 8)
}, numeric(n2)))
rownames(bg_counts2) <- bg_genes

sex_counts2 <- matrix(0, nrow = length(sex_genes), ncol = n2, dimnames = list(sex_genes, sample_id2))
is_female2 <- sex2 == "Female"
sex_counts2["XIST", is_female2]  <- rnbinom(sum(is_female2),  mu = 900 * lib_factor2[is_female2],  size = 10)
sex_counts2["XIST", !is_female2] <- rnbinom(sum(!is_female2), mu = 8   * lib_factor2[!is_female2],  size = 10)
for (g in setdiff(sex_genes, "XIST")) {
  sex_counts2[g, is_female2]  <- rnbinom(sum(is_female2),  mu = 5   * lib_factor2[is_female2],  size = 10)
  sex_counts2[g, !is_female2] <- rnbinom(sum(!is_female2), mu = 400 * lib_factor2[!is_female2], size = 10)
}

expected_count2 <- rbind(sex_counts2, bg_counts2)[gene_symbols, ]
storage.mode(expected_count2) <- "double"
colnames(expected_count2) <- sample_id2

tpm2 <- sweep(expected_count2, 2, colSums(expected_count2), FUN = function(x, s) x / s * 1e6)
fpkm2 <- tpm2

se2 <- SummarizedExperiment(
  assays = list(expected_count = expected_count2, TPM = tpm2, FPKM = fpkm2),
  rowData = row_data,
  colData = DataFrame(qc_cols2)
)

# Deliberately NO tissue_type/collection_site/cohort columns here -- this is the
# whole point of the second scenario.
rnaAnnot2 <- data.frame(
  sample_id = sample_id2,
  patient_id = patient_id2,
  Adequacy_group = adequacy_group,
  batch = batch2,
  rin = rin2,
  dv200 = dv200_2,
  stringsAsFactors = FALSE
)

pairs2 <- combn(sample_id2, 2)
n_pairs2 <- ncol(pairs2)
hets_a2 <- round(runif(n_pairs2, 2000, 5000))
hom_alts_a2 <- round(hets_a2 * runif(n_pairs2, 0.3, 0.7))
hets_b2 <- round(runif(n_pairs2, 2000, 5000))
hom_alts_b2 <- round(hets_b2 * runif(n_pairs2, 0.3, 0.7))

somalier_pairs_out2 <- data.frame(
  "#sample_a" = pairs2[1, ],
  sample_b = pairs2[2, ],
  relatedness = round(rbeta(n_pairs2, 1.5, 8), 3),
  hets_a = hets_a2,
  hom_alts_a = hom_alts_a2,
  hets_b = hets_b2,
  hom_alts_b = hom_alts_b2,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

genotypes_out2 <- data.frame(subject = sample_id2, stringsAsFactors = FALSE)
for (locus in hla_loci) {
  pool <- allele_pool(locus)
  genotypes_out2[[paste0(locus, "_1")]] <- sample(pool, n2, replace = TRUE)
  genotypes_out2[[paste0(locus, "_2")]] <- sample(pool, n2, replace = TRUE)
}

out_dir2 <- file.path(repo_dir, "example_data", "data2_style")
dir.create(out_dir2, recursive = TRUE, showWarnings = FALSE)

se_path2 <- file.path(out_dir2, "Aggregated_Gene_Expression.rds")
rnaAnnot_path2 <- file.path(out_dir2, "metadata.tsv")
somalier_pairs_path2 <- file.path(out_dir2, "somalier.pairs.tsv")
genotypes_path2 <- file.path(out_dir2, "genotypes.tsv")

saveRDS(se2, se_path2)
write.table(rnaAnnot2, rnaAnnot_path2, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(somalier_pairs_out2, somalier_pairs_path2, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(genotypes_out2, genotypes_path2, sep = "\t", row.names = FALSE, quote = FALSE)

se2_in <- readRDS(se_path2)
rnaAnnot2_in <- read.table(rnaAnnot_path2, sep = "\t", header = TRUE)
somalier_pairs2_in <- read.delim(somalier_pairs_path2, header = TRUE, stringsAsFactors = FALSE)
genotypes2_in <- read.delim(genotypes_path2)

output_dir2 <- file.path(out_dir2, "output")
dir.create(output_dir2, recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  file.path(repo_dir, "qc_report.Rmd"),
  envir = new.env(parent = globalenv()),
  params = list(
    se = se2_in,
    rnaAnnot = rnaAnnot2_in,
    somalier_pairs = somalier_pairs2_in,
    genotypes = genotypes2_in,
    # No tissue_type/collection_site/cohort entries -- they simply don't exist in
    # this cohort's metadata, and column_map silently skips anything not present.
    column_map = list(
      sample_id = "sample_id", patient_id = "patient_id", batch_id = "batch",
      rin = "rin", dv200 = "dv200"
    ),
    tissueType = "Data2-style Synthetic Example QC",
    tissue = NULL,
    group_vars = c("Adequacy_group"),
    output_se = file.path(output_dir2, "annotated_Gene_Expression.rds"),
    output_se_per = file.path(output_dir2, "all_samples_Gene_Expression.rds"),
    qcFile = file.path(output_dir2, "qc_flag_summary.tsv"),
    showSession = TRUE
  ),
  output_file = "example_QCReport_data2style.html",
  output_dir = output_dir2,
  intermediates_dir = output_dir2
)

message("Done. Data2-style report written to: ", file.path(output_dir2, "example_QCReport_data2style.html"))

source(file.path(repo_dir, "generate_somalier_network.R"))
generate_somalier_network(
  somalier_pairs = somalier_pairs2_in,
  output_html = file.path(output_dir2, "somalier_network_graph.html"),
  cutoff = 0.6
)
