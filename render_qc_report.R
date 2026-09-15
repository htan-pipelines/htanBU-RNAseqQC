#!/usr/bin/env Rscript
# ==============================================================================
# render_qc_report.R
#
# Example/template driver for htanBU-RNAseqQC (qc_report.Rmd).
#
# Renders qc_report.Rmd once per tissue_type found in your metadata, after both
# stages of the bulk-rna-seq-pipeline have been run:
#   1. RNA_seq_pipeline.wdl  (per sample)
#   2. Aggregation/aggregate.wdl (across the cohort)
#      https://github.com/htan-pipelines/bulk-rna-seq-pipeline
#
# Copy this file, edit the CONFIG block below to point at your own files, and run:
#   Rscript render_qc_report.R
# ==============================================================================

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(rmarkdown)
})

# ------------------------------------------------------------------------------
# CONFIG -- edit everything in this block for your own project
# ------------------------------------------------------------------------------
config <- list(
  # Path to the aggregated SummarizedExperiment .rds from Aggregation/aggregate.wdl
  se_path = "path/to/your/Aggregated_Gene_Expression.rds",

  # Path to your external per-sample metadata table (tab-separated)
  rnaAnnot_path = "path/to/your/metadata.tsv",

  # Path to somalier pairwise relatedness output from Aggregation/aggregate.wdl
  somalier_pairs_path = "path/to/your/somalier.pairs.tsv",

  # Path to arcasHLA genotype calls from Aggregation/aggregate.wdl
  genotypes_path = "path/to/your/genotypes.tsv",

  # Map YOUR metadata column names onto the canonical names the report expects.
  # Must match params$column_map in qc_report.Rmd. tissue_type/collection_site/cohort/
  # batch_id are all optional -- if your cohort doesn't have one of these concepts,
  # just omit that entry (or leave it NULL).
  column_map = list(
    sample_id       = "sample_id",         # must match colnames(se) after mapping
    patient_id      = "patient_id",
    tissue_type     = "tissue_type",
    collection_site = "collection_site",
    cohort          = "cohort",
    batch_id        = "batch_id",
    rin             = "rin",
    dv200           = "dv200"
  ),

  # Free-text label used in the report title/header
  report_title = "Bulk RNA-seq QC",

  # Which column to split reports on. Set to NULL for a single-tissue cohort (or any
  # cohort you don't want split by tissue) -- one report covers everything.
  tissue_col = "tissue_type",

  # Extra categorical breakdown variable(s) (as they appear in rnaAnnot after
  # column_map) for the general-stats plots, PCA shape aesthetic, and somalier
  # heatmap annotation -- see params$group_vars in qc_report.Rmd. Default reproduces
  # the original collection_site + cohort breakdowns; point this at your own
  # variable(s) of interest instead (e.g. c("Adequacy_group")) if your cohort doesn't
  # have those.
  group_vars = c("collection_site", "cohort"),

  # QC flagging cutoffs -- see qc_report.Rmd YAML header for definitions.
  # Leave as-is to use the report's built-in defaults, or override per project.
  qc_cutoffs = list(
    tin_median = 50, rin = 5, threeprime_bias = 0.5, genes_detected = 5000,
    exon_cv = 1.0, rrna_rate = 0.01, dv200 = 50, heterozygosity_mean = 0.5,
    somalier_relatedness = 0.6, hla_match = 0.6
  ),

  # Optional: spreadsheet highlighting a priority sample subset (NULL to skip)
  priority_list_file = NULL,

  # Output locations
  output_dir = file.path("QC_output", format(Sys.Date(), "%Y%m%d")),
  save_annotated_se = TRUE   # write the full annotated SE (output_se param) once
)

# ------------------------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------------------------
se <- readRDS(config$se_path)
rnaAnnot <- read.table(config$rnaAnnot_path, sep = "\t", header = TRUE)
somalier_pairs <- read.delim(config$somalier_pairs_path, header = TRUE, stringsAsFactors = FALSE)
genotypes <- read.delim(config$genotypes_path)



dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Optionally save the full annotated SE once, up front
# ------------------------------------------------------------------------------
output_se <- if (isTRUE(config$save_annotated_se)) {
  file.path(config$output_dir, paste0(format(Sys.Date(), "%Y%m%d"), "_annotated_Gene_Expression.rds"))
} else NULL

# ------------------------------------------------------------------------------
# Derive the tissue values to report on directly from the metadata, instead of
# hardcoding a fixed list -- this works for any cohort/tissue set. NULL config$
# tissue_col means "don't split by tissue" -- render_one() is then called once,
# covering the whole cohort.
# ------------------------------------------------------------------------------
tissue_values <- if (!is.null(config$tissue_col)) {
  sort(unique(as.character(rnaAnnot[[config$tissue_col]])))
} else {
  NULL
}

render_one <- function(tissue_value = NULL) {
  tag <- if (!is.null(tissue_value)) gsub("[^A-Za-z0-9]+", "", tissue_value) else "all_samples"
  out_subdir <- file.path(config$output_dir, tag)
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

  message("Rendering QC report for ", if (!is.null(tissue_value)) paste0("tissue='", tissue_value, "'") else "the whole cohort", " ...")

  rmarkdown::render(
    "qc_report.Rmd",
    params = list(
      se = se,
      rnaAnnot = rnaAnnot,
      somalier_pairs = somalier_pairs,
      genotypes = genotypes,
      column_map = config$column_map,
      tissueType = config$report_title,
      tissue = tissue_value,
      group_vars = config$group_vars,
      qc_cutoffs = config$qc_cutoffs,
      priority_list_file = config$priority_list_file,
      output_se = output_se,
      output_se_per = file.path(out_subdir, paste0(format(Sys.Date(), "%Y%m%d"), "_", tag, "_Gene_Expression.rds")),
      qcFile = file.path(out_subdir, paste0(format(Sys.Date(), "%Y%m%d"), "_", tag, "_QC_FlagSummary.tsv")),
      showSession = TRUE
    ),
    output_file = paste0(format(Sys.Date(), "%Y%m%d"), "_", tag, "_QCReport.html"),
    output_dir = out_subdir,
    intermediates_dir = out_subdir
  )
}

if (is.null(tissue_values)) {
  render_one(NULL)
} else {
  for (tv in tissue_values) {
    render_one(tv)
  }
}

# ------------------------------------------------------------------------------
# Also generate the interactive somalier relatedness network as a standalone
# HTML file (kept out of the Rmd report itself -- see generate_somalier_network.R)
# ------------------------------------------------------------------------------
source("generate_somalier_network.R")
generate_somalier_network(
  somalier_pairs = somalier_pairs,
  output_html = file.path(config$output_dir, "somalier_network_graph.html"),
  cutoff = config$qc_cutoffs$somalier_relatedness
)

message("Done. Reports written under: ", normalizePath(config$output_dir))
