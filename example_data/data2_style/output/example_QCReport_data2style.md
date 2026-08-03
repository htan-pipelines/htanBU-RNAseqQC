---
title: "htanBU-RNAseqQC: Bulk RNA-seq QC Report"
date: "2026-08-03"
output:
 html_document:
     toc: true
     toc_float: true
     number_sections: true
     theme: cosmo
     code_folding: hide
     keep_md: true
params:
  # ---- Required inputs (no sensible default; render() will error early if missing) ----
  se: NULL                 # SummarizedExperiment produced by aggregate.wdl (combine_se task)
  rnaAnnot: NULL            # data.frame/tibble of external per-sample metadata (e.g. Terra data table export)
  somalier_pairs: NULL      # data.frame from Aggregation/somalier_final.wdl pairwise output (somalier.pairs.tsv)
  genotypes: NULL           # data.frame from Aggregation/arcasHLA_merge.wdl (genotypes.tsv)

  # ---- Column-name mapping ----
  # This report standardizes on the canonical column names below (right-hand side).
  # Point each entry at whatever column holds that information in YOUR rnaAnnot file.
  # Only `sample_id` is strictly required to match colnames(se); `tissue_type`,
  # `collection_site`, and `cohort` are all fully optional -- if your cohort doesn't
  # have one of these concepts (e.g. a single-tissue study), just omit that entry (or
  # leave it NULL) and the corresponding sections/plots are skipped automatically.
  # `group_vars` (below) doesn't need a column_map entry at all -- any column present
  # in rnaAnnot after this mapping, including ones with no canonical slot here (e.g. a
  # study-specific "Adequacy_group"), can be named directly.
  column_map: !r list(
    sample_id       = "entity.pcga2_biospecimen_id",
    patient_id      = "patient_ID",
    tissue_type     = "tissue_type",
    collection_site = "collection_site",
    cohort          = "cohort",
    batch_id        = "sequencing_batch_id",
    rin             = "rin",
    dv200           = "dv200"
    )

  # ---- Report labeling / subsetting ----
  tissueType: "Bulk RNA-seq QC"   # free-text title used in report headers
  tissue: NULL                     # value of `tissue_type` to subset on for section 2 (NULL = use all samples,
                                    # or if there's no `tissue_type` column at all)

  # ---- Extra grouping/breakdown variables ----
  # Column name(s) (post column_map, as they appear in rnaAnnot) to use as additional
  # categorical breakdown dimensions in the general-stats plots, PCA shape aesthetic,
  # and somalier heatmap annotation -- wherever `collection_site`/`cohort` show up
  # below. Default reproduces the original collection_site + cohort breakdowns. For a
  # cohort without those concepts, point this at whatever variable(s) of interest your
  # study actually has, e.g. `c("Adequacy_group")`.
  group_vars: !r c("collection_site", "cohort")

  # ---- QC flagging cutoffs (see "QC cutoffs used in this report" for definitions) ----
  qc_cutoffs: !r list(
    tin_median            = 50,
    rin                   = 5,
    threeprime_bias       = 0.5,
    genes_detected        = 5000,
    exon_cv                = 1.0,
    rrna_rate              = 0.01,
    dv200                  = 50,
    heterozygosity_mean    = 0.5,
    somalier_relatedness   = 0.6,
    hla_match               = 0.6
    )

  # ---- Optional inputs ----
  priority_list_file: NULL   # optional .xlsx with a "Priority" column + ID columns to highlight a sample subset;
                              # the "Priority samples" heatmap section is skipped when this is NULL
  priority_id_cols: !r c("RNA biospecimen ID", "Parent biospecimen ID")
  priority_column: "Priority"
  priority_value: "High"

  # ---- Outputs ----
  output_se: NULL          # path to save the full annotated SE (skipped if NULL)
  output_se_per: NULL      # path to save the tissue-subset SE (skipped if NULL)
  qcFile: "qc_flag_summary.tsv"   # path to write the per-sample QC flag table (TSV)
  showSession: TRUE
---

# Data2-style Synthetic Example QC

This report performs quality control (QC) on bulk RNA-seq data produced by the
[bulk-rna-seq-pipeline](https://github.com/htan-pipelines/bulk-rna-seq-pipeline)
(`RNA_seq_pipeline.wdl` run per sample, aggregated with `Aggregation/aggregate.wdl`).
It reads the aggregated `SummarizedExperiment`, layers on external sample metadata,
somalier relatedness, and HLA genotype calls, and walks through: general sample
statistics, sample-swap/contamination checks, per-metric QC scatterplots against
configurable cutoffs, expression-based outlier detection (correlation + PCA),
QC-metric clustering, sex-check from XIST/Y-linked gene expression, and a combined
per-sample QC flag summary.

## Required inputs

This report is meant to run **after** you have completed both stages of the
bulk-rna-seq-pipeline:

1. **`RNA_seq_pipeline.wdl`** run once per sample (produces the per-sample gene/isoform
   `SummarizedExperiment`, FastQC, STAR, RNA-SeQC2, RSeQC TIN, samtools, somalier
   extraction, and arcasHLA genotyping outputs — see pipeline diagram).
2. **`Aggregation/aggregate.wdl`** run once across all samples in the cohort (produces
   the combined `SummarizedExperiment`, `somalier.pairs.tsv`, and `genotypes.tsv`).

You need to supply the following four objects/files as `params` when rendering:

| Param | What it is | Produced by | Example |
|---|---|---|---|
| `se` | Aggregated `SummarizedExperiment` (assays: `expected_count`, `TPM`, `FPKM`; `colData` has per-sample QC columns prefixed `fastqc_`, `STAR_`, `rnaseqc_`, `samtools_`, `TIN_`, plus `Somalier`/`Genotypes` summaries) | `Aggregation/aggregate.wdl` -> `combine_se.wdl` | `PCGA2_Gene_Expression.rds`, loaded with `readRDS()` |
| `rnaAnnot` | Per-sample external metadata table (sample ID, patient/participant ID, tissue type, collection site, cohort, batch, RIN, DV200, etc.) | Your own sample-tracking sheet / Terra data table export | `read.table("metadata.tsv", sep = "\t", header = TRUE)` |
| `somalier_pairs` | Pairwise sample relatedness | `Aggregation/aggregate.wdl` -> `somalier_final.wdl` | `read.delim("somalier.pairs.tsv")` |
| `genotypes` | Per-sample HLA genotype calls | `Aggregation/aggregate.wdl` -> `arcasHLA_merge.wdl` | `read.delim("genotypes.tsv")` |

The `column_map` param tells the report which columns in `rnaAnnot` correspond to the
sample ID, patient ID, tissue type, collection site, cohort, and batch — update it to
match your own metadata sheet's column names; you do not need to rename anything in your
source file. `sample_id` (post-mapping) must match `colnames(se)` exactly.

See `render_qc_report.R` in this repo for a complete, runnable example, and
`generate_somalier_network.R` for the interactive sample-relatedness network
(rendered to a standalone `somalier_network_graph.html`, not embedded in this report).

## Setup
Load required libraries, read in parameters, validate required inputs, apply the
column-name mapping, and set up display options for the report.


```r
# NOTE on dependencies: the original PCGA02 script loaded 33 packages, of which
# roughly half (DESeq2, sva, patchwork, gridExtra, lubridate, stringr, stringi,
# viridis, fgsea, ggprism, rstatix, glue, Hmisc, mixtools, plotmm, fields,
# ggraph, jsonlite) were never actually called anywhere in the report code --
# most likely carried over from a shared lab template. They've been dropped
# here to keep the dependency footprint (and install time) down for a public
# QC pipeline; the list below is only what qc_report.Rmd actually uses.
# `Matrix`, `knitr`, `methods`, and `utils` are used via `::` and don't need
# an explicit library() call (Matrix must still be installed).
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(edgeR)
  library(RColorBrewer)
  library(circlize)
  library(ggplot2)
  library(pheatmap)
  library(grid)
  library(readxl)
  library(ggpubr)
  library(ComplexHeatmap)
  library(ggrepel)
  library(plotly)
  library(tidyverse)
  library(igraph)
  library(kableExtra)
})

# ---- Read in params ----
se <- params$se
rnaAnnot <- params$rnaAnnot
somalier_pairs <- params$somalier_pairs
genotypes <- params$genotypes
column_map <- params$column_map
tissueType <- params$tissueType
tissue <- params$tissue
group_vars <- params$group_vars
qc_cutoffs <- params$qc_cutoffs
priority_list_file <- params$priority_list_file
priority_id_cols <- params$priority_id_cols
priority_column <- params$priority_column
priority_value <- params$priority_value
output_se <- params$output_se
output_se_per <- params$output_se_per
qcFile <- params$qcFile
showSession <- params$showSession

# ---- Validate required inputs up front, fail fast with a clear message ----
required_params <- c("se", "rnaAnnot", "somalier_pairs", "genotypes")
missing_params <- required_params[vapply(mget(required_params), is.null, logical(1))]
if (length(missing_params) > 0) {
  stop(
    "htanBU-RNAseqQC: missing required param(s): ", paste(missing_params, collapse = ", "),
    ". See the 'Required inputs' section of this report / README.md for what each one is."
  )
}
if (!methods::is(se, "SummarizedExperiment")) {
  stop("htanBU-RNAseqQC: `se` must be a SummarizedExperiment object (see combine_se.wdl output).")
}

# Small helper used throughout the report to build a discrete color palette that always
# has exactly as many colors as there are factor levels (avoids silently recycled /
# missing colors when a cohort has more or fewer groups than the original PCGA dataset).
mk_disc_cols <- function(x, palette = "Set1") {
  lev <- sort(unique(as.character(x)))
  n <- length(lev)
  base <- RColorBrewer::brewer.pal(min(9, max(3, n)), palette)
  cols <- if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
  setNames(cols, lev)
}

# plotly::ggplotly() can throw ("subscript out of bounds" in gg2list's plot-margin
# unit conversion, among other internal errors) on some ggplot objects depending on
# plotly/ggplot2 version and the specific plot's scale/legend configuration -- a
# real, observed failure mode on real cohort data that's been hard to pin to one
# root cause. Since these plots are also perfectly viewable as static ggplots, don't
# let a plotly conversion failure take down the whole report: fall back to the
# static plot (with a warning) instead of erroring out.
safe_ggplotly <- function(p, ...) {
  tryCatch(
    plotly::ggplotly(p, ...),
    error = function(e) {
      warning("ggplotly() failed (", conditionMessage(e), "); showing a static plot instead.")
      p
    }
  )
}

# Set up Rmarkdown display options
dev <- c("png", "pdf")
knitr::opts_chunk$set(
  echo = TRUE,
  warning = FALSE,
  message = FALSE,
  cache = FALSE,
  cache.lazy = FALSE,
  cache.comments = FALSE,
  fig.align = "center",
  fig.keep = "all",
  dev = dev
)
```

## Clean up metadata

Apply `column_map` so the rest of the report can refer to canonical column names
(`sample_id`, `patient_id`, `tissue_type`, `collection_site`, `cohort`, `batch_id`,
`RIN`, `DV200`) regardless of what they're called in your source metadata file.


```r
# Rename whichever columns in rnaAnnot are configured in column_map to their
# canonical names. Missing/optional fields (e.g. no `cohort` column) are simply skipped.
for (canonical_name in names(column_map)) {
  source_name <- column_map[[canonical_name]]
  if (!is.null(source_name) && source_name %in% colnames(rnaAnnot)) {
    colnames(rnaAnnot)[colnames(rnaAnnot) == source_name] <- canonical_name
  }
}
colnames(rnaAnnot)[colnames(rnaAnnot) == "rin"] <- "RIN"
colnames(rnaAnnot)[colnames(rnaAnnot) == "dv200"] <- "DV200"

if (!"sample_id" %in% colnames(rnaAnnot)) {
  stop("htanBU-RNAseqQC: after applying `column_map`, rnaAnnot has no `sample_id` column. ",
       "Update params$column_map$sample_id to the column in your metadata file that matches colnames(se).")
}

mismatched_ids <- rnaAnnot$sample_id[!(rnaAnnot$sample_id %in% colnames(se))]
if (length(mismatched_ids) > 0) {
  warning("htanBU-RNAseqQC: ", length(mismatched_ids),
          " sample_id(s) in rnaAnnot were not found in colnames(se): ",
          paste(utils::head(mismatched_ids, 10), collapse = ", "),
          if (length(mismatched_ids) > 10) ", ..." else "")
}

# `tissue_type` is fully optional: cohorts with a single tissue (or none tracked at
# all) simply skip every tissue_type-specific section/subsetting step below.
has_tissue_type <- "tissue_type" %in% colnames(rnaAnnot)
if (!has_tissue_type) tissue <- NULL

# Drop any group_vars entries that don't actually exist in rnaAnnot (e.g. a typo, or
# a column_map mismatch) rather than erroring deep inside some later plot.
bad_group_vars <- group_vars[!group_vars %in% colnames(rnaAnnot)]
if (length(bad_group_vars) > 0) {
  warning("htanBU-RNAseqQC: ", length(bad_group_vars),
          " params$group_vars not found in rnaAnnot (after column_map): ",
          paste(bad_group_vars, collapse = ", "), " -- dropping.")
}
group_vars <- setdiff(group_vars, bad_group_vars)

# Re-order rnaAnnot to match the sample order in the SummarizedExperiment
rnaAnnot <- rnaAnnot[match(colnames(se), rnaAnnot$sample_id), ]
rownames(rnaAnnot) <- gsub(" ", "", rnaAnnot$sample_id)
```


## General Statistics Plots {.tabset .tabset-fade}


```r
# Always runs (regardless of which of tissue_type/group_vars are present) since the
# per-group_var panels below all depend on it. Built from whatever combination of
# patient_id + tissue_type (if present) + group_vars actually exists in rnaAnnot.
core_cols <- c("patient_id", if (has_tissue_type) "tissue_type", group_vars)
summary.samples <- rnaAnnot %>%
  dplyr::select(sample_id, dplyr::all_of(core_cols)) %>%
  dplyr::group_by(dplyr::across(dplyr::all_of(core_cols))) %>%
  dplyr::count(patient_id)
```


```r
# results='asis' + a conditional cat() heading (matching the group_var_panels
# pattern below) so the tab is skipped ENTIRELY -- not just left empty -- when
# there's no tissue_type column. A plain `eval=` on this chunk would still leave
# behind a static "### Tissue type" markdown heading with nothing under it, since
# that heading lives outside the chunk.
if (has_tissue_type) {
  cat("\n\n### Tissue type\n\n")

  colors <- mk_disc_cols(rnaAnnot$tissue_type, "Set1")

  p <- ggplot(rnaAnnot, aes(tissue_type, fill = tissue_type)) +
    geom_bar(position = position_fill()) +
    geom_text(aes(label = ..count..),
              stat = "count",
              colour = "white",
              position = position_fill(vjust = 0.3),size=9) +
    geom_text(aes(label = paste0(round(..count../sum(..count..)*100,2), "%")),
              stat = "count",
              colour = "white",
              position = position_fill(vjust = 0.5),size=9) +
    scale_fill_manual(values = colors) +
    coord_flip() +
    guides(fill = guide_legend(title = "Tissue type"))+
    theme_bw() %+replace%
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            text = element_text(size = 19, family = "serif"),
            axis.ticks = element_blank(),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color = "black"),
            legend.position = "top",
            legend.text = element_text(size = 15),
            legend.key.size = unit(1, "char")
      )
  print(p)
  cat("\n")
}
```


```r
if (has_tissue_type) {
  cat("\n\n### Tissue type per participant\n\n")

  colors <- mk_disc_cols(rnaAnnot$tissue_type, "Set1")

  p <- ggplot(summary.samples, aes(x=patient_id, y=n))+
    geom_bar(aes(fill=tissue_type),stat="identity",position="stack")+
    labs(y="Total tissue count", x = "Participant")+
    scale_y_continuous(breaks = round(seq(min(summary.samples$n),15, by = 1),1))+
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(title = "Tissue type"))+
  theme_bw() %+replace%
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            text = element_text(size = 19, family = "serif"),
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle=-90, hjust=.1),
            axis.text.y = element_text(color = "black"),
            legend.position = "top",
            legend.text = element_text(size = 20),
            legend.key.size = unit(1, "char")
      )
  print(p)
  cat("\n")
}
```


```r
# One tab per params$group_vars entry (default c("collection_site", "cohort"),
# reproducing the original two fixed "per site"/"per cohort" panels exactly). When
# tissue_type is present, each panel breaks that variable down by tissue_type, same as
# before; otherwise it falls back to a simple frequency count of the variable alone.
for (v in group_vars) {
  cat("\n\n### ", v, "\n\n")

  colors <- mk_disc_cols(rnaAnnot[[v]], "Set1")
  x_var <- if (has_tissue_type) "tissue_type" else v

  p <- ggplot(summary.samples, aes(x = .data[[x_var]], y = n, fill = .data[[v]])) +
    geom_bar(stat = "identity", position = "stack") +
    labs(y = "Total count", x = if (has_tissue_type) NULL else v) +
    scale_fill_manual(values = colors) +
    guides(fill = guide_legend(title = v)) +
    theme_bw() %+replace%
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            text = element_text(size = 19, family = "serif"),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_text(color = "black"),
            legend.position = "top",
            legend.text = element_text(size = 20),
            legend.key.size = unit(1, "char")
      )

  print(p)
  cat("\n")
}
```



###  Adequacy_group 

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/group_var_panels-1.png" style="display: block; margin: auto;" />


```r
# All cutoffs come from params$qc_cutoffs (see YAML header) so they can be
# overridden per cohort/tissue without editing this file.
# NOTE: this chunk used to live just above "QC Metric Overview" (after tissue
# subsetting), but the Sample Swap/Contamination Check section above (HLA genotype
# check specifically) reads `hla_cutoff` -- moved up here, right after setup, so
# every cutoff is available before it's first used anywhere in the report.
tin_cutoff <- qc_cutoffs$tin_median
rin_cutoff <- qc_cutoffs$rin
threeprime_bias_cutoff <- qc_cutoffs$threeprime_bias
gene_detected_cutoff <- qc_cutoffs$genes_detected
exon_cv_cutoff <- qc_cutoffs$exon_cv
rrna_rate_cutoff <- qc_cutoffs$rrna_rate
dv200_cutoff <- qc_cutoffs$dv200
hetero_cutoff <- qc_cutoffs$heterozygosity_mean
somalier_cutoff <- qc_cutoffs$somalier_relatedness
hla_cutoff <- qc_cutoffs$hla_match
```

## Sample Swap/Contamination Check {.tabset .tabset-fade}
Create pairwise matrix of somalier sample pairs output, calculate heterozygosity rate, and library complexity to help identify poor quality samples.


```r
sample_id <- rnaAnnot$sample_id
participant_id <- rnaAnnot$patient_id
# Initialize matrix
mat.sf <- matrix(0, nrow = length(sample_id), ncol = length(sample_id), dimnames = list(sample_id, sample_id))

# Iterate over sample pairs and populate matrix
for (i in seq_along(sample_id)) {
  relatedness_values <- somalier_pairs %>% 
    filter(sample_b == sample_id[i]) %>% 
    pull(relatedness)  # Extracts relatedness as a vector
  
  if (length(relatedness_values) > 0) {  # Check if non-empty
    mat.sf[i, 1:length(relatedness_values)] <- relatedness_values
  }
}

# Make matrix symmetric
mat.sf.sym <- as.matrix(Matrix::forceSymmetric(mat.sf, uplo="L"))

# Assign row and column names
rownames(mat.sf.sym) <- sample_id
colnames(mat.sf.sym) <- sample_id


####### Compute heterozygosity rates
somalier_pairs <- somalier_pairs %>%
  mutate(
    heterozygosity_rate_a = round(hets_a / (hets_a + hom_alts_a), 3),
    heterozygosity_rate_b = round(hets_b / (hets_b + hom_alts_b), 3)
  )

# Initialize matrices
het_data <- matrix(0, nrow = length(sample_id), ncol = length(sample_id), dimnames = list(sample_id, sample_id))
het_mat <- matrix(0, nrow = length(sample_id), ncol = length(sample_id), dimnames = list(sample_id, sample_id))

# Populate matrices
for (i in seq_along(sample_id)) {
  # Extract heterozygosity rate A
  het_values_a <- somalier_pairs %>%
    filter(X.sample_a == sample_id[i]) %>%
    pull(heterozygosity_rate_a)  # Extract as vector

  # Extract heterozygosity rate B
  het_values_b <- somalier_pairs %>%
    filter(sample_b == sample_id[i]) %>%
    pull(heterozygosity_rate_b)

  if (length(het_values_a) > 0) het_data[i, 1:length(het_values_a)] <- het_values_a
  if (length(het_values_b) > 0) het_mat[i, 1:length(het_values_b)] <- het_values_b
}

# Combine heterozygosity data
het_comb <- het_data[, 1]
het_comb[length(het_comb)] <- het_mat[nrow(het_mat), 1]

# Add heterozygosity mean to RNA annotation
rnaAnnot$heterozygosity_mean <- het_comb

####### Compute Library complexity 
# #Gene count data is the expected counts stored in the summarized experiment
#Create a DGE object
data.dge <- c()
data.dge <- DGEList(counts = assays(se)$expected_count)

## rank genes with most reads to fewest reads and take the top 50 genes
indx.50 =apply(data.dge$counts, 2, function(x) order(x, decreasing = TRUE))[1:50]
top50genes <- data.dge$counts[indx.50,]
lib.compelx.50 <- diag(apply(data.dge$counts, 2, function(x) (apply(top50genes, 2, function(x) sum(x)))/ sum(x))) # output values diag

## rank genes with most reads to fewest reads and take the top 500 genes
indx.500 =apply(data.dge$counts, 2, function(x) order(x, decreasing = TRUE))[1:500]
top500genes <- data.dge$counts[indx.500,]
lib.compelx.500 <- diag(apply(data.dge$counts, 2, function(x) (apply(top500genes, 2, function(x) sum(x)))/ sum(x))) # output values diag

## rank genes with most reads to fewest reads and take the top 1000 genes
indx.1000 =apply(data.dge$counts, 2, function(x) order(x, decreasing = TRUE))[1:1000]
top1000genes <- data.dge$counts[indx.1000,]
lib.compelx.1000 <- diag(apply(data.dge$counts, 2, function(x) (apply(top1000genes, 2, function(x) sum(x)))/ sum(x))) # output values diag

rnaAnnot$library_complexity_50 <- lib.compelx.50
rnaAnnot$library_complexity_500 <- lib.compelx.500
rnaAnnot$library_complexity_1000 <- lib.compelx.1000
#add rna annotation to the summarized experiment
colData(se) <- cbind(colData(se),rnaAnnot[,colnames(rnaAnnot)[!colnames(rnaAnnot) %in% colnames(colData(se))]])
```


### Filtering out non-numeric QC metrics from SE object 


```r
select_qc_metrics <- function(se_object) {
  #Remove fastqc and other QC metrics 
  list_to_remove = c("fastqc","Somalier","Genotypes","rnaseqc_Sample",colnames(rnaAnnot)[!colnames(rnaAnnot) %in% c("RIN","DV200","heterozygosity_mean")])
  #Empty vector to store index of qc metrics 
  qc_to_exclude = c()
  for (i in list_to_remove){
    names_qc <- grep(paste0("^",i),colnames(colData(se_object)),ignore.case = TRUE)
    qc_to_exclude <- c(qc_to_exclude,names_qc)
    #remove_t <- grep(i,colnames(colData(se)),ignore.case = TRUE)
  }
  #Remove fastqc and other QC metrics and keep the rest 
  qc.metrics <- colData(se_object)[-qc_to_exclude]
  #Remove % sign for some rows 
  qc.metrics@listData <- lapply(qc.metrics@listData, FUN = function(x){gsub("%","",x)})
  # convert character type to numeric type 
  qc.metrics@listData <- lapply(qc.metrics@listData, FUN = function(y){if(is.character(y)) as.numeric(y) else y})
  return(qc.metrics)
}
## All samples 
qc.metrics <- select_qc_metrics(se_object = se)
```

### Somalier relatedness Heatmap semi-supervised 
 

```r
##### heatmap---- start
mat.sf.sym <- apply(mat.sf.sym, 1, function(x) ifelse(x > 1, 1, ifelse(x < -1, -1, x)))
#mat.sf.sym <- apply(mat.sf.sym, 1, function(x) ifelse(x > 0.4, 0.4, ifelse(x < -0.4, 0, x)))

# mk_disc_cols() is defined once in the setup chunk and reused here. Discrete
# annotation variables are tissue_type (if present) + params$group_vars -- dynamic so
# this works regardless of which of those a given cohort actually has.
disc_vars <- c(if (has_tissue_type) "tissue_type", group_vars)
colorList <- setNames(
  lapply(disc_vars, function(v) mk_disc_cols(rnaAnnot[[v]], "Set1")),
  disc_vars
)

cont_fun <- function(v) circlize::colorRamp2(
  quantile(v, c(0.05, 0.5, 0.95), na.rm = TRUE),
  c("navy", "white", "firebrick")
)
# continuous annotation columns
colorList <- c(colorList, list(
  star_uniquely_mapped_per             = cont_fun(qc.metrics$STAR_uniquely_mapped_percent),
  rnaseqc_median_transcript_coverage   = cont_fun(qc.metrics$rnaseqc_Median.of.Avg.Transcript.Coverage),
  rnaseqc_rRNA_rate                    = cont_fun(qc.metrics$rnaseqc_rRNA.Rate),
  rnaseqc_high_quality_reads           = cont_fun(qc.metrics$rnaseqc_High.Quality.Reads),
  samtools_properly_paired             = cont_fun(qc.metrics$samtools_reads_properly_paired),
  ranseqc_low_quality_reads            = cont_fun(qc.metrics$rnaseqc_Low.Quality.Reads)
))

# NOTE: as.data.frame(lapply(character(0), ...)) is a 0-row/0-col data.frame, not an
# N-row one -- cbind()-ing that against N-row data.frames below would error, so build
# a proper N-row-but-0-column placeholder when there are no discrete vars at all.
disc_df <- if (length(disc_vars) > 0) {
  setNames(as.data.frame(lapply(disc_vars, function(v) rnaAnnot[[v]])), disc_vars)
} else {
  data.frame(row.names = seq_len(nrow(rnaAnnot)))
}

ha = HeatmapAnnotation(
  df = cbind(
    data.frame(participant = rnaAnnot$patient_id),
    disc_df,
    data.frame(
      star_uniquely_mapped_per = qc.metrics$STAR_uniquely_mapped_percent,
      rnaseqc_median_transcript_coverage = qc.metrics$rnaseqc_Median.of.Avg.Transcript.Coverage,
      rnaseqc_rRNA_rate = qc.metrics$rnaseqc_rRNA.Rate,
      rnaseqc_high_quality_reads = qc.metrics$rnaseqc_High.Quality.Reads,
      samtools_properly_paired = qc.metrics$samtools_reads_properly_paired,
      ranseqc_low_quality_reads = qc.metrics$rnaseqc_Low.Quality.Reads
    )
  ),
  annotation_height = unit(4, "mm"),annotation_name_side = "left", col = colorList
)
dend1 = cluster_within_group(mat.sf.sym, rnaAnnot$patient_id)
Heatmap(mat.sf.sym, name = "relatedness", top_annotation = ha,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(at = c(-1, 0, 1)),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),
        cluster_columns = dend1, cluster_rows = TRUE, column_names_centered = T
        #column_names_rot = 45,
        #row_names_max_width = unit(10, "cm"),
        #column_names_max_height = unit(9,"cm"),
        #rowAnnotation(hetero_rate = anno_barplot(het_comb, width = unit(4, "cm")))
) +
  Heatmap(het_comb, name = "heterozygosity rate", row_names_gp = gpar(fontsize = 5),
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6),
                                                                    height = unit(2, "cm"),
                                                                    width = unit(2,"cm"))), width = unit(7, "mm")) +
  rowAnnotation("heterozygosity rate" = anno_barplot(het_comb, width = unit(4, "cm"),gp = gpar(fill = 5)))
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/somalier_semi_supervised_heatmap-1.png" style="display: block; margin: auto;" />

### Somalier relatedness Heatmap unsupervised 


```r
# pdf("/rprojectnb/pulmseq/SFStudy/qc/somalier_heatmaps/heatmap_somalier_relatedness_QC_cluster.pdf", width = 24, height = 12)
Heatmap(mat.sf.sym, name = "relatedness", top_annotation = ha,
        show_row_names = FALSE, show_column_names = TRUE,
        heatmap_legend_param = list(at = c(-1, 0, 1)),
        row_names_gp = gpar(fontsize = 5),
        column_names_gp = gpar(fontsize = 5),
        #cluster_columns = TRUE, cluster_rows = TRUE, column_names_centered = T
        cluster_columns = TRUE, cluster_rows = TRUE, column_names_centered = T
        #column_names_rot = 45,
        #row_names_max_width = unit(10, "cm"),
        #column_names_max_height = unit(9,"cm"),
        #rowAnnotation(hetero_rate = anno_barplot(het_comb, width = unit(4, "cm")))
) +
  Heatmap(het_comb, name = "heterozygosity rate", row_names_gp = gpar(fontsize = 5),
          top_annotation = HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 2:6),
                                                                    height = unit(2, "cm"),
                                                                    width = unit(2,"cm"))), width = unit(7, "mm")) +
  rowAnnotation("heterozygosity rate" = anno_barplot(het_comb, width = unit(4, "cm"),gp = gpar(fill = 5)))
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/somalier_unsupervised_heatmap-1.png" style="display: block; margin: auto;" />

```r
# dev.off()
```

### Somalier relatedness network graph

The interactive sample-relatedness network (colored by patient ID and by TIN) is
generated as a **standalone HTML file** rather than embedded in this report — run
`generate_somalier_network.R` (in this repo) against your `somalier.pairs.tsv` to
produce `somalier_network_graph.html`. It uses the same `somalier_pairs` relatedness cutoff
(`params$qc_cutoffs$somalier_relatedness`, default 0.6) as the flagging logic below,
and lets you interactively search/filter samples by patient — something a static plot
embedded in this report can't do.

<!--
The `hla_folder` (per-sample sample.genotype.json / sample.genes.json from arcasHLA
genotype) based cross-checking approach was dropped from this report: it duplicated
the section below (which uses the already-aggregated genotypes.tsv from
arcasHLA_merge.wdl) and required a per-sample file layout most users won't have handy.
If you do have the per-sample JSON files and want read-count-weighted HLA comparisons,
see the git history of PCGA02_BulkRNA_QC_v2.Rmd in the PCGA project repo for the
original (commented-out) implementation.
-->

### HLA Genotype Quality Assessment (use genotypes.tsv)

```r
genotype.pairs <- function(x) {
  # Generate all pairwise sample combinations
  pairs <- combn(rownames(x), 2, simplify = FALSE)

  # Preallocate result as a data frame
  comp <- data.frame(
    sample_a = character(length(pairs)),
    sample_b = character(length(pairs)),
    percentOFcomparison = numeric(length(pairs)),
    stringsAsFactors = FALSE
  )

  # Compute pairwise comparisons efficiently
  for (i in seq_along(pairs)) {
    row_a <- pairs[[i]][1]
    row_b <- pairs[[i]][2]

    a_vec <- x[row_a, ]
    b_vec <- x[row_b, ]

    # Compute similarity as a fraction of total HLA fields being compared
    # (e.g. arcasHLA's default 7-gene x 2-allele genotype call = 14 columns;
    # computed from ncol(x) rather than hardcoded so this works for any
    # HLA-typing panel width)
    shared_features <- sum(a_vec %in% b_vec)
    comp$sample_a[i] <- row_a
    comp$sample_b[i] <- row_b
    comp$percentOFcomparison[i] <- round(shared_features / ncol(x), 2)
  }

  return(comp)
}

# Preprocess the genotypes param (already read in by render_qc_report.R / your caller;
# NOTE: this used to re-read a hardcoded "TerraOutput/genotypes.tsv" path here, silently
# ignoring params$genotypes -- fixed so the report always uses what was passed in)
genotypes <- genotypes %>%
  arrange(subject) %>%
  column_to_rownames("subject")  # Convert "subject" column to rownames
geno.all.pairs <- genotype.pairs(genotypes)

# Create an adjacency matrix
samples <- unique(c(geno.all.pairs$sample_a, geno.all.pairs$sample_b))
n <- length(samples)
adj_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(samples, samples))

# Iterate over each row in geno.all.pairs
for (i in 1:nrow(geno.all.pairs)) {
  sample_a <- geno.all.pairs$sample_a[i]  # First sample
  sample_b <- geno.all.pairs$sample_b[i]  # Second sample
  matchPercent <- geno.all.pairs$percentOFcomparison[i]  # Relatedness value

  # Only add edge if relatedness is >= the configured HLA match cutoff
  if (matchPercent >= hla_cutoff) {
    idx1 <- which(rownames(adj_matrix) == sample_a)
    idx2 <- which(rownames(adj_matrix) == sample_b)

    # Ensure both samples exist in merged_qc before adding edges
    if (length(idx1) > 0 & length(idx2) > 0) {
      adj_matrix[idx1, idx2] <- 1  # Add edge
      adj_matrix[idx2, idx1] <- 1  # Ensure symmetry
    }
  }
}

# Create an igraph object from the adjacency matrix
network <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = NULL, diag = FALSE)


patient_ids <- sapply(strsplit(samples,"_"), function(x) paste(x[1:2], collapse = "_"))
patient_ids_uni <- unique(patient_ids)
colors <- hcl.colors(length(patient_ids_uni), palette = "Set3")
patient_color_map <- setNames(colors, patient_ids_uni)
node_colors_patient <- patient_color_map[patient_ids]

# Layout with Fruchterman-Reingold to avoid label overlap
# layout <- layout_with_fr(network)
layout <- layout_with_fr(network, niter = 5000, grid = "nogrid")


plot(network,
     layout = layout, # Use Fruchterman-Reingold layout to separate nodes
     vertex.size = 3, # Adjust vertex size
     vertex.label.cex = 0.8,  # Adjust label size
     vertex.label.color = "black",  # Color of labels
     vertex.label.family = "sans",  # Font for labels
     vertex.color = node_colors_patient,
     main = "Sample Network Based on Confident HLA Genotype, Cutoff = 0.6,Colored by Patient ID",
     vertex.label.dist = 0.5,  # Increase label distance to avoid overlap
     edge.color = "black",  # Color of edges
     edge.width = 1)      # Width of edges
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/HLA_genotype_check_alt-1.png" style="display: block; margin: auto;" />

```r
legend("topright",
       legend = patient_ids_uni,
       fill = colors,
       title = "Patient",
       cex = 0.8)
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/HLA_genotype_check_alt-2.png" style="display: block; margin: auto;" />


# QC: (Data2-style Synthetic Example QC) Analysis
## QC Metric Overview {.tabset .tabset-fade}


```r
# tissue_type/tissue are optional -- with neither, se.subset is just the whole cohort.
se.subset <- if (has_tissue_type && !is.null(tissue)) {
  se[, which(colData(se)$tissue_type == tissue)]
} else {
  se
}
qcMetrics.subset <- select_qc_metrics(se_object = se.subset)
```


### TIN vs RIN

```r
df <- as.data.frame(qcMetrics.subset)
df$sample_id <- rownames(df)

flag <- subset(df, TIN_median < tin_cutoff & RIN < rin_cutoff)

ggscatter(df, x="TIN_median", y="RIN",
          add="reg.line", conf.int=TRUE,
          cor.coef=TRUE, cor.method="pearson",
          xlab="Transcript Integrity Number (TIN)",
          ylab="RNA Integrity Number (RIN)",
          color="black", repel=FALSE) +
  geom_text_repel(
    data = flag,
    aes(TIN_median, RIN, label = sample_id),
    color = "red"
  )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/TIN_vs_RIN-1.png" style="display: block; margin: auto;" />
Correlation between Transcript integrity number (TIN) and RNA integrity number (RIN).


### TIN vs Median 3 bias 


```r
df <- as.data.frame(qcMetrics.subset)
df$sample_id <- rownames(df)

flag <- subset(df, rnaseqc_Median.3..bias > threeprime_bias_cutoff)

ggscatter(
  df,
  x = "TIN_median",
  y = "rnaseqc_Median.3..bias",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "pearson",
  xlab = "Transcript Integrity Number (TIN)",
  ylab = "Median 3' Bias",
  repel = TRUE
) +
  geom_text_repel(
    data = flag,
    aes(x = TIN_median, y = rnaseqc_Median.3..bias, label = sample_id),
    color = "red"
  )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/TIN_vs_Median_3_bias-1.png" style="display: block; margin: auto;" />

Correlation between Transcript integrity number (TIN) and Median 3 bias.


### Genes Detected vs Median Exon CV 


```r
df <- as.data.frame(qcMetrics.subset)
df$sample_id <- rownames(df)

flag <- subset(
  df,
  rnaseqc_Genes.Detected < gene_detected_cutoff |
  rnaseqc_Median.Exon.CV > exon_cv_cutoff
)

ggscatter(
  df,
  x = "rnaseqc_Genes.Detected",
  y = "rnaseqc_Median.Exon.CV",
  color = "rnaseqc_Duplicate.Rate.of.Mapped",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "pearson",
  xlab = "Genes Detected",
  ylab = "Median Exon CV",
  repel = TRUE
) +
  geom_text_repel(
    data = flag,
    aes(x = rnaseqc_Genes.Detected, y = rnaseqc_Median.Exon.CV, label = sample_id),
    color = "red"
  )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/Genes_detected_vs_Median_Exon_CV-1.png" style="display: block; margin: auto;" />

Plot shows the coefficient of variation (CV) of exon coverage facilitates identification of lower-quality samples with higher variability in coverage, higher duplication rates and fewer genes detected.

### RIN vs rRNA rate 


```r
df <- as.data.frame(qcMetrics.subset)
df$sample_id <- rownames(df)

flag <- subset(df, rnaseqc_rRNA.Rate > rrna_rate_cutoff)

ggscatter(
  df,
  x = "RIN",
  y = "rnaseqc_rRNA.Rate",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "pearson",
  xlab = "RNA Integrity Number (RIN)",
  ylab = "rRNA Rate",
  repel = TRUE
) +
  geom_text_repel(
    data = flag,
    aes(x = RIN, y = rnaseqc_rRNA.Rate, label = sample_id),
    color = "red"
  )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/RIN_vs_rRNA_rate-1.png" style="display: block; margin: auto;" />

Correlation between RNA integrity number (RIN) and percentage of of ribosomal RNA.


### TIN vs DV200


```r
df <- as.data.frame(qcMetrics.subset)
df$sample_id <- rownames(df)

flag <- subset(df, DV200 < dv200_cutoff)

ggscatter(
  df,
  x = "TIN_median",
  y = "DV200",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "pearson",
  xlab = "Transcript Integrity Number (TIN)",
  ylab = "DV200",
  repel = TRUE
) +
  geom_text_repel(
    data = flag,
    aes(x = TIN_median, y = DV200, label = sample_id),
    color = "red"
  )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/TIN_vs_DV200-1.png" style="display: block; margin: auto;" />

Correlation between Transcript integrity number (TIN) and DV200.

### RIN vs DV200

```r
df <- as.data.frame(qcMetrics.subset)
df$sample_id <- rownames(df)

# label low-quality samples (adjust or remove if not needed)
flag <- subset(df, DV200 < dv200_cutoff)

ggscatter(
  df,
  x = "RIN",
  y = "DV200",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "pearson",
  xlab = "RNA integrity number (RIN)",
  ylab = "DV200",
  repel = TRUE
) +
  geom_text_repel(
    data = flag,
    aes(x = RIN, y = DV200, label = sample_id),
    color = "red"
  )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/RIN_vs_DV200-1.png" style="display: block; margin: auto;" />

Correlation between RNA integrity number (RIN) and DV200.

### Somalier Heterozygosity vs DV200

```r
df <- as.data.frame(qcMetrics.subset)
df$sample_id <- rownames(df)

flag <- subset(df, heterozygosity_mean > hetero_cutoff)

ggscatter(
  df,
  x = "heterozygosity_mean",
  y = "DV200",
  add = "reg.line", conf.int = TRUE,
  cor.coef = TRUE, cor.method = "pearson",
  xlab = "Somalier Heterozygosity Mean",
  ylab = "DV200",
  repel = TRUE
) +
  geom_text_repel(
    data = flag,
    aes(x = heterozygosity_mean, y = DV200, label = sample_id),
    color = "red"
  )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/Somalier_heterozygosity_vs_DV200-1.png" style="display: block; margin: auto;" />




```r
# Normalize qc metrics for heatmap 
normalized_qc_metrics <- function(metric_data) {
  # Z-score normalize the quality data for plotting via heatmap
  normalized.data <- t(scale(t(t(as.matrix(metric_data))),center=T,scale=T))
  # Eliminate any rows that have NA after normalization
  normalized.data<-normalized.data[which(is.na(normalized.data[,1])==FALSE),]
  return(normalized.data)
}

normalizedQC <- normalized_qc_metrics(qcMetrics.subset)
```

### QC Metrics Heatmap 

```r
#QC of interest 
qc_to_keep <- c("rnaseqc_Total.Reads",
                "rnaseqc_Mapped.Reads", 
                "rnaseqc_Mapped.Unique.Reads",
                "rnaseqc_Unique.Rate.of.Mapped",
                "rnaseqc_High.Quality.Reads", 
                "rnaseqc_rRNA.Reads",
                "rnaseqc_rRNA.Rate",
                "rnaseqc_Median.3..bias",
                "rnaseqc_Genes.Detected",
                "rnaseqc_Median.Exon.CV",
                "rnaseqc_Duplicate.Rate.of.Mapped",
                "RIN", 
                "TIN_median", 
                "STAR_uniquely_mapped_percent", 
                "STAR_num_splices",
                "STAR_unmapped_tooshort_percent",
                "STAR_avg_mapped_read_length",
                "STAR_uniquely_mapped",
                "STAR_total_reads", 
                "STAR_multimapped_multiple",
                "STAR_multimapped_multiple_percent",
                "DV200")

#Subset QC of interest 
normalizedQC.sub <- normalizedQC[rownames(normalizedQC) %in% qc_to_keep,]

normalizedQC.sub <- t(apply(normalizedQC.sub, 1, function(x) ifelse(x > 6, 6, ifelse(x < -6, -6, x))))


# Annotate by tissue_type if present, else fall back to the first group_var, else
# skip annotation entirely (pheatmap's own default for "no annotation" is NA, not NULL).
annot_var <- if (has_tissue_type) "tissue_type" else if (length(group_vars) > 0) group_vars[1] else NA
if (!is.na(annot_var)) {
  mat_col <- data.frame(sample = as.factor(colData(se.subset)[[annot_var]]))
  rownames(mat_col) <- colnames(se.subset)
} else {
  mat_col <- NA
}


#Plot a heatmap of QC metrics across all samples
heatmap.qc <- pheatmap(
mat=normalizedQC.sub,
col=colorRampPalette(c("blue","white","red"), space="rgb")(255),
annotation_col = mat_col,
clustering_method = "ward.D2",
fontsize_row = 10,
legend = TRUE,
legend_breaks = c(-6,0,6),
border_color = FALSE,
show_colnames = TRUE)

heatmap.qc
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/qc_metrics_heatmap-1.png" style="display: block; margin: auto;" />


## QC Based on Expression Count {.tabset .tabset-fade}


```r
# ------------------------------------
# 1. Expression Count Preprocessing
# ------------------------------------

# Create a DGE object from expression counts in SummarizedExperiment
data.dge <- c()
data.dge <- DGEList(counts = assays(se.subset)$expected_count)

# Normalize library sizes using TMM method
data.dge <- calcNormFactors(data.dge, method = "TMM")

# Compute log-CPM values
data.cpm <- cpm(data.dge, log = TRUE)

# ------------------------------------
# 2. Filter Out Lowly Expressed Genes
# ------------------------------------

# Identify lowly expressed genes based on variance and expression levels
test1 <- apply(data.cpm, 1, IQR)  # Interquartile range (variance)
test2 <- apply(data.cpm, 1, sum)  # Sum of CPM values
low.genes <- union(which(test1 == 0), which(test2 <= 1))  # Zero variance or CPM < 1
gene.ind <- setdiff(1:nrow(data.cpm), low.genes)  # Retain expressed genes

# Recompute normalization on filtered genes
data.dge2 <- data.dge[gene.ind, , keep.lib.sizes = FALSE]
data.dge2 <- calcNormFactors(data.dge2, method = "TMM")
data.cpm2 <- cpm(data.dge2, log = TRUE)

# Z-score normalization for heatmap
# The scale() function will always scale by column, only (you can get it to scale by row by doing t(scale(t(x)))); so, each column in the data is scaled separately.
z.data <- t(scale(t(data.cpm2), center = TRUE, scale = TRUE))
```


### Sample-Sample correlaton heatmap with RIN values


```r
#Plots sample correlation heatmaps
sampleMatrix <- cor(data.cpm2)
gene.cor.tally <- matrix(ncol = nrow(sampleMatrix), data = 0)
colnames(gene.cor.tally) <- rownames(sampleMatrix)

gene.cor.tally2 <- matrix(
  ncol = nrow(sampleMatrix),
  nrow = ncol(sampleMatrix),
  data = 0)
rownames(gene.cor.tally2) <- colnames(sampleMatrix)
colnames(gene.cor.tally2) <- rownames(sampleMatrix)

isNum <- 0
for (i in 1:ncol(sampleMatrix)){
  qc <- sampleMatrix[,i]
  isNum <- isNum + 1
    mean.qc <- mean(qc)
    sd.qc <- sd(qc)
    # Identify samples ±2 SD from mean correlation
    outlier_samples <- which(qc > mean.qc + 2 * sd.qc | qc < mean.qc - 2 * sd.qc)
    gene.cor.tally[,outlier_samples] <- gene.cor.tally[,outlier_samples] + 1
    
    gene.cor.tally2[i,outlier_samples] <- gene.cor.tally2[i, outlier_samples] + 1
}
# samples that are +/- 2sd away from the mean sample to sample correlation 
colnames(gene.cor.tally)[gene.cor.tally > 10]
```

```
## character(0)
```

Outliers are defined as +/- 2sd away from the mean sample to sample correlation of all genes

This plot shows hierarchically clustering of pearson correlation of log cpm transformed expression data. 

```r
rownames(sampleMatrix) <- paste(colnames(se.subset),colData(se.subset)$RIN, sep="-")
colnames(sampleMatrix) <- colnames(se.subset) #paste(colnames(se),colData(se)$tissue_type, sep="-")
colors <- colorRampPalette(c("blue", "white", "red"))(255)

# Create annotation dataframe for RIN values
annotation_rin <- data.frame(RIN = qcMetrics.subset$RIN)
rownames(annotation_rin) <- rownames(qcMetrics.subset)

# Heatmap with RIN annotation
pheatmap(sampleMatrix,
         fontsize = 8,
         color = colors,
         annotation_col = annotation_rin,
         clustering_method = "ward.D2",
         border_color = FALSE,
         show_colnames = TRUE,
         main = "Sample Correlation Heatmap (RIN Median Annotated)")
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/Sample_sample_correlation_heatmap_RIN-1.png" style="display: block; margin: auto;" />

### Sample-Sample correlaton heatmap with TIN values


```r
rownames(sampleMatrix) <- paste(colnames(se.subset),round(colData(se.subset)$TIN_median,3), sep="-")
colnames(sampleMatrix) <- colnames(se.subset)#paste(colnames(se),colData(se)$tissue_type, sep="-")

# Create annotation dataframe for RIN values
annotation_tin <- data.frame(TIN = qcMetrics.subset$TIN_median)
rownames(annotation_tin) <- rownames(qcMetrics.subset)

# Heatmap with RIN annotation
pheatmap(sampleMatrix,
         fontsize = 8,
         color = colors,
         annotation_col = annotation_tin,
         clustering_method = "ward.D2",
         border_color = FALSE,
         show_colnames = TRUE,
         main = "Sample Correlation Heatmap (TIN Median Annotated)")
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/Sample_sample_correlation_heatmap_TIN-1.png" style="display: block; margin: auto;" />

Similar pattern is observed as above, 
### Outlier sample names based on PCA


```r
#Conduct a PCA analysis and plot PC1 versus PC2
pca <- prcomp(z.data,scale=F, center=F)
mean.pc1 <- mean(pca$rotation[,1])
sd.pc1 <- sd(pca$rotation[,1])
mean.pc2 <- mean(pca$rotation[,2])
sd.pc2 <- sd(pca$rotation[,2])
outliers <- c(which(pca$rotation[,1]<(mean.pc1-sd.pc1*2)),
              which(pca$rotation[,1]>(mean.pc1+sd.pc1*2)),
              which(pca$rotation[,2]<(mean.pc2-sd.pc2*2)),
              which(pca$rotation[,2]>(mean.pc2+sd.pc2*2)))
out.color <- rep("no", length(colData(se.subset)$sample_id))
out.color[outliers] <- "yes"
perc.var <- (pca$sdev)^2/sum(pca$sdev^2)
percentVar <- round(100 * perc.var)
# sample name of outliers 
colData(se.subset)$sample_id[which(out.color == "yes")]
```

```
## character(0)
```
Outlier samples defined +/- 2sd away from the mean PC1 ad mean PC2

### PCA plot colored by TIN


```r
# Shape by the first group_var (default collection_site), if any are configured.
shape_var <- if (length(group_vars) > 0) group_vars[1] else NULL
sample_id <- rownames(pca$rotation)
plot2 <- ggplot(as.data.frame(pca$rotation[,c(1,2)]),
                aes(PC1,PC2,color=colData(se.subset)$TIN_median,label = sample_id)) +
  (if (!is.null(shape_var)) geom_point(size=3,aes(shape=colData(se.subset)[[shape_var]])) else geom_point(size=3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(shape = shape_var)
safe_ggplotly(plot2)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-28e95e8fcf3ee2dd04cc" style="width:1056px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-28e95e8fcf3ee2dd04cc">{"x":{"data":[{"x":[-0.34132555669682618,0.19930779475461491,-0.31978116431201509,0.26440980228685385,-0.21468381307301776,0.26409980068254812],"y":[0.002143706201169712,-0.54228177375048936,-0.013105187639462749,0.14215137559030044,0.4837524741811246,-0.16606817341186514],"text":["PC1: -0.3413256<br />PC2:  0.002143706<br />colData(se.subset)$TIN_median: 60.2<br />sample_id: Data2Cohort_P01_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1:  0.1993078<br />PC2: -0.542281774<br />colData(se.subset)$TIN_median: 45.2<br />sample_id: Data2Cohort_P02_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1: -0.3197812<br />PC2: -0.013105188<br />colData(se.subset)$TIN_median: 70.7<br />sample_id: Data2Cohort_P05_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1:  0.2644098<br />PC2:  0.142151376<br />colData(se.subset)$TIN_median: 45.7<br />sample_id: Data2Cohort_P06_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1: -0.2146838<br />PC2:  0.483752474<br />colData(se.subset)$TIN_median: 80.7<br />sample_id: Data2Cohort_P09_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1:  0.2640998<br />PC2: -0.166068173<br />colData(se.subset)$TIN_median: 53.9<br />sample_id: Data2Cohort_P10_S1<br />colData(se.subset)[[shape_var]]: Adequate"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(41,86,125,1)","rgba(19,43,67,1)","rgba(57,119,169,1)","rgba(20,44,69,1)","rgba(73,152,213,1)","rgba(31,67,100,1)"],"opacity":1,"size":11.338582677165356,"symbol":"circle","line":{"width":1.8897637795275593,"color":["rgba(41,86,125,1)","rgba(19,43,67,1)","rgba(57,119,169,1)","rgba(20,44,69,1)","rgba(73,152,213,1)","rgba(31,67,100,1)"]}},"hoveron":"points","name":"Adequate","legendgroup":"Adequate","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.33147941528946495,-0.19479135720787935,-0.29838270612673901],"y":[-0.16292302179947601,0.36308795005655825,-0.41814741950596496],"text":["PC1: -0.3314794<br />PC2: -0.162923022<br />colData(se.subset)$TIN_median: 88.2<br />sample_id: Data2Cohort_P03_S1<br />colData(se.subset)[[shape_var]]: Borderline","PC1: -0.1947914<br />PC2:  0.363087950<br />colData(se.subset)$TIN_median: 83.7<br />sample_id: Data2Cohort_P07_S1<br />colData(se.subset)[[shape_var]]: Borderline","PC1: -0.2983827<br />PC2: -0.418147420<br />colData(se.subset)$TIN_median: 51.4<br />sample_id: Data2Cohort_P11_S1<br />colData(se.subset)[[shape_var]]: Borderline"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(86,177,247,1)","rgba(78,162,226,1)","rgba(28,60,90,1)"],"opacity":1,"size":11.338582677165356,"symbol":"triangle-up","line":{"width":1.8897637795275593,"color":["rgba(86,177,247,1)","rgba(78,162,226,1)","rgba(28,60,90,1)"]}},"hoveron":"points","name":"Borderline","legendgroup":"Borderline","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33961684010506166,0.2763920102131891,0.35661776466367789],"y":[-0.043233468011760666,0.063338233143606607,0.29128530494625321],"text":["PC1:  0.3396168<br />PC2: -0.043233468<br />colData(se.subset)$TIN_median: 54.7<br />sample_id: Data2Cohort_P04_S1<br />colData(se.subset)[[shape_var]]: Intermediate","PC1:  0.2763920<br />PC2:  0.063338233<br />colData(se.subset)$TIN_median: 56.5<br />sample_id: Data2Cohort_P08_S1<br />colData(se.subset)[[shape_var]]: Intermediate","PC1:  0.3566178<br />PC2:  0.291285305<br />colData(se.subset)$TIN_median: 66.9<br />sample_id: Data2Cohort_P12_S1<br />colData(se.subset)[[shape_var]]: Intermediate"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(33,70,103,1)","rgba(35,75,110,1)","rgba(51,107,153,1)"],"opacity":1,"size":11.338582677165356,"symbol":"square","line":{"width":1.8897637795275593,"color":["rgba(33,70,103,1)","rgba(35,75,110,1)","rgba(51,107,153,1)"]}},"hoveron":"points","name":"Intermediate","legendgroup":"Intermediate","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.20000000000000001],"y":[-0.30000000000000004],"name":"b2d8307290a7932aa211a8139a2fbd43","type":"scatter","mode":"markers","opacity":0,"hoverinfo":"skip","showlegend":false,"marker":{"color":[0,1],"colorscale":[[0,"#132B43"],[0.003344481605351184,"#132B44"],[0.0066889632107023679,"#132C44"],[0.010033444816053552,"#142C45"],[0.013377926421404736,"#142D45"],[0.016722408026755918,"#142D46"],[0.020066889632107104,"#142D46"],[0.023411371237458123,"#142E47"],[0.026755852842809305,"#152E47"],[0.030100334448160491,"#152F48"],[0.033444816053511676,"#152F48"],[0.036789297658862859,"#152F49"],[0.040133779264214041,"#153049"],[0.043478260869565223,"#16304A"],[0.046822742474916405,"#16304A"],[0.050167224080267594,"#16314B"],[0.053511705685618777,"#16314B"],[0.056856187290969959,"#16324C"],[0.060200668896320982,"#17324D"],[0.063545150501672157,"#17324D"],[0.066889632107023353,"#17334E"],[0.070234113712374535,"#17334E"],[0.073578595317725717,"#17344F"],[0.0769230769230769,"#18344F"],[0.080267558528428082,"#183450"],[0.083612040133779264,"#183550"],[0.086956521739130446,"#183551"],[0.090301003344481628,"#183651"],[0.093645484949832811,"#193652"],[0.096989966555184007,"#193652"],[0.10033444816053519,"#193753"],[0.10367892976588637,"#193754"],[0.10702341137123739,"#193854"],[0.11036789297658857,"#1A3855"],[0.11371237458193975,"#1A3955"],[0.11705685618729093,"#1A3956"],[0.12040133779264212,"#1A3956"],[0.12374581939799331,"#1A3A57"],[0.12709030100334448,"#1B3A57"],[0.13043478260869568,"#1B3B58"],[0.13377926421404684,"#1B3B59"],[0.13712374581939804,"#1B3B59"],[0.14046822742474907,"#1C3C5A"],[0.14381270903010024,"#1C3C5A"],[0.14715719063545143,"#1C3D5B"],[0.1505016722408026,"#1C3D5B"],[0.1538461538461538,"#1C3D5C"],[0.15719063545150497,"#1D3E5C"],[0.16053511705685616,"#1D3E5D"],[0.16387959866220736,"#1D3F5D"],[0.16722408026755853,"#1D3F5E"],[0.17056856187290972,"#1D3F5F"],[0.17391304347826089,"#1E405F"],[0.17725752508361209,"#1E4060"],[0.18060200668896326,"#1E4160"],[0.18394648829431445,"#1E4161"],[0.18729096989966562,"#1E4261"],[0.19063545150501665,"#1F4262"],[0.19397993311036785,"#1F4263"],[0.19732441471571902,"#1F4363"],[0.20066889632107021,"#1F4364"],[0.20401337792642138,"#1F4464"],[0.20735785953177258,"#204465"],[0.21070234113712374,"#204465"],[0.21404682274247494,"#204566"],[0.21739130434782614,"#204566"],[0.22073578595317714,"#214667"],[0.22408026755852833,"#214668"],[0.2274247491638795,"#214768"],[0.2307692307692307,"#214769"],[0.23411371237458187,"#214769"],[0.23745819397993306,"#22486A"],[0.24080267558528423,"#22486A"],[0.24414715719063543,"#22496B"],[0.24749163879598662,"#22496C"],[0.25083612040133779,"#224A6C"],[0.25418060200668896,"#234A6D"],[0.25752508361204018,"#234A6D"],[0.26086956521739135,"#234B6E"],[0.26421404682274252,"#234B6E"],[0.26755852842809369,"#244C6F"],[0.27090301003344475,"#244C70"],[0.27424749163879591,"#244C70"],[0.27759197324414708,"#244D71"],[0.28093645484949831,"#244D71"],[0.28428093645484948,"#254E72"],[0.28762541806020064,"#254E72"],[0.29096989966555181,"#254F73"],[0.29431438127090304,"#254F74"],[0.2976588628762542,"#254F74"],[0.30100334448160521,"#265075"],[0.30434782608695643,"#265075"],[0.3076923076923076,"#265176"],[0.31103678929765877,"#265176"],[0.31438127090300994,"#275277"],[0.31772575250836116,"#275278"],[0.32107023411371233,"#275278"],[0.3244147157190635,"#275379"],[0.32775919732441472,"#275379"],[0.33110367892976589,"#28547A"],[0.33444816053511706,"#28547B"],[0.33779264214046822,"#28557B"],[0.34113712374581945,"#28557C"],[0.34448160535117062,"#28567C"],[0.34782608695652178,"#29567D"],[0.35117056856187295,"#29567D"],[0.35451505016722401,"#29577E"],[0.35785953177257518,"#29577F"],[0.36120401337792635,"#2A587F"],[0.36454849498327757,"#2A5880"],[0.36789297658862874,"#2A5980"],[0.37123745819397991,"#2A5981"],[0.37458193979933108,"#2A5982"],[0.37792642140468213,"#2B5A82"],[0.3812709030100333,"#2B5A83"],[0.38461538461538447,"#2B5B83"],[0.38795986622073569,"#2B5B84"],[0.39130434782608686,"#2C5C85"],[0.39464882943143803,"#2C5C85"],[0.3979933110367892,"#2C5D86"],[0.40133779264214042,"#2C5D86"],[0.40468227424749159,"#2C5D87"],[0.40802675585284276,"#2D5E87"],[0.41137123745819398,"#2D5E88"],[0.41471571906354515,"#2D5F89"],[0.41806020066889632,"#2D5F89"],[0.42140468227424749,"#2E608A"],[0.42474916387959871,"#2E608A"],[0.42809364548494988,"#2E618B"],[0.43143812709030105,"#2E618C"],[0.43478260869565227,"#2E618C"],[0.43812709030100344,"#2F628D"],[0.44147157190635461,"#2F628D"],[0.44481605351170578,"#2F638E"],[0.448160535117057,"#2F638F"],[0.45150501672240817,"#30648F"],[0.45484949832775934,"#306490"],[0.45819397993311051,"#306590"],[0.4615384615384614,"#306591"],[0.46488294314381257,"#306592"],[0.46822742474916373,"#316692"],[0.47157190635451496,"#316693"],[0.47491638795986613,"#316793"],[0.47826086956521729,"#316794"],[0.48160535117056846,"#326895"],[0.48494983277591969,"#326895"],[0.48829431438127086,"#326996"],[0.49163879598662202,"#326996"],[0.49498327759197325,"#326997"],[0.49832775919732442,"#336A98"],[0.50167224080267558,"#336A98"],[0.50501672240802675,"#336B99"],[0.50836120401337792,"#336B99"],[0.51170568561872909,"#346C9A"],[0.51505016722408037,"#346C9B"],[0.51839464882943154,"#346D9B"],[0.52173913043478237,"#346D9C"],[0.52508361204013354,"#346E9D"],[0.52842809364548471,"#356E9D"],[0.53177257525083588,"#356E9E"],[0.53511705685618705,"#356F9E"],[0.53846153846153832,"#356F9F"],[0.54180602006688949,"#3670A0"],[0.54515050167224066,"#3670A0"],[0.54849498327759183,"#3671A1"],[0.551839464882943,"#3671A1"],[0.55518394648829417,"#3772A2"],[0.55852842809364533,"#3772A3"],[0.56187290969899661,"#3773A3"],[0.56521739130434778,"#3773A4"],[0.56856187290969895,"#3773A4"],[0.57190635451505012,"#3874A5"],[0.57525083612040129,"#3874A6"],[0.57859531772575246,"#3875A6"],[0.58193979933110362,"#3875A7"],[0.5852842809364549,"#3976A8"],[0.58862876254180607,"#3976A8"],[0.59197324414715724,"#3977A9"],[0.59531772575250841,"#3977A9"],[0.59866220735785958,"#3978AA"],[0.60200668896321075,"#3A78AB"],[0.60535117056856191,"#3A79AB"],[0.60869565217391319,"#3A79AC"],[0.61204013377926436,"#3A79AC"],[0.61538461538461553,"#3B7AAD"],[0.6187290969899667,"#3B7AAE"],[0.62207357859531753,"#3B7BAE"],[0.6254180602006687,"#3B7BAF"],[0.62876254180601987,"#3C7CB0"],[0.63210702341137115,"#3C7CB0"],[0.63545150501672232,"#3C7DB1"],[0.63879598662207349,"#3C7DB1"],[0.64214046822742465,"#3C7EB2"],[0.64548494983277582,"#3D7EB3"],[0.64882943143812699,"#3D7FB3"],[0.65217391304347816,"#3D7FB4"],[0.65551839464882944,"#3D7FB5"],[0.65886287625418061,"#3E80B5"],[0.66220735785953178,"#3E80B6"],[0.66555183946488294,"#3E81B6"],[0.66889632107023411,"#3E81B7"],[0.67224080267558528,"#3F82B8"],[0.67558528428093645,"#3F82B8"],[0.67892976588628773,"#3F83B9"],[0.6822742474916389,"#3F83BA"],[0.68561872909698973,"#4084BA"],[0.6889632107023409,"#4084BB"],[0.69230769230769207,"#4085BB"],[0.69565217391304324,"#4085BC"],[0.69899665551839441,"#4086BD"],[0.70234113712374568,"#4186BD"],[0.70568561872909685,"#4186BE"],[0.70903010033444802,"#4187BF"],[0.71237458193979919,"#4187BF"],[0.71571906354515036,"#4288C0"],[0.71906354515050153,"#4288C1"],[0.72240802675585269,"#4289C1"],[0.72575250836120386,"#4289C2"],[0.72909698996655514,"#438AC2"],[0.73244147157190631,"#438AC3"],[0.73578595317725748,"#438BC4"],[0.73913043478260865,"#438BC4"],[0.74247491638795982,"#438CC5"],[0.74581939799331098,"#448CC6"],[0.74916387959866215,"#448DC6"],[0.75250836120401343,"#448DC7"],[0.7558528428093646,"#448EC8"],[0.75919732441471577,"#458EC8"],[0.76254180602006694,"#458FC9"],[0.76588628762541811,"#458FC9"],[0.76923076923076927,"#458FCA"],[0.77257525083612044,"#4690CB"],[0.77591973244147172,"#4690CB"],[0.77926421404682289,"#4691CC"],[0.78260869565217406,"#4691CD"],[0.78595317725752523,"#4792CD"],[0.78929765886287639,"#4792CE"],[0.79264214046822756,"#4793CF"],[0.7959866220735784,"#4793CF"],[0.79933110367892968,"#4894D0"],[0.80267558528428085,"#4894D0"],[0.80602006688963201,"#4895D1"],[0.80936454849498318,"#4895D2"],[0.81270903010033435,"#4896D2"],[0.81605351170568552,"#4996D3"],[0.81939799331103669,"#4997D4"],[0.82274247491638797,"#4997D4"],[0.82608695652173914,"#4998D5"],[0.8294314381270903,"#4A98D6"],[0.83277591973244147,"#4A99D6"],[0.83612040133779231,"#4A99D7"],[0.83946488294314348,"#4A9AD8"],[0.84280936454849464,"#4B9AD8"],[0.84615384615384592,"#4B9BD9"],[0.84949832775919709,"#4B9BDA"],[0.85284280936454826,"#4B9BDA"],[0.85618729096989943,"#4C9CDB"],[0.8595317725752506,"#4C9CDB"],[0.86287625418060176,"#4C9DDC"],[0.86622073578595293,"#4C9DDD"],[0.86956521739130421,"#4D9EDD"],[0.87290969899665538,"#4D9EDE"],[0.87625418060200655,"#4D9FDF"],[0.87959866220735772,"#4D9FDF"],[0.88294314381270889,"#4DA0E0"],[0.88628762541806005,"#4EA0E1"],[0.88963210702341122,"#4EA1E1"],[0.8929765886287625,"#4EA1E2"],[0.89632107023411367,"#4EA2E3"],[0.89966555183946484,"#4FA2E3"],[0.90301003344481601,"#4FA3E4"],[0.90635451505016718,"#4FA3E5"],[0.90969899665551834,"#4FA4E5"],[0.91304347826086951,"#50A4E6"],[0.91638795986622068,"#50A5E7"],[0.91973244147157196,"#50A5E7"],[0.92307692307692313,"#50A6E8"],[0.9264214046822743,"#51A6E8"],[0.92976588628762546,"#51A7E9"],[0.93311036789297663,"#51A7EA"],[0.9364548494983278,"#51A8EA"],[0.93979933110367897,"#52A8EB"],[0.94314381270903025,"#52A9EC"],[0.94648829431438142,"#52A9EC"],[0.94983277591973259,"#52AAED"],[0.95317725752508375,"#53AAEE"],[0.95652173913043492,"#53ABEE"],[0.95986622073578576,"#53ABEF"],[0.96321070234113693,"#53ACF0"],[0.96655518394648821,"#54ACF0"],[0.96989966555183937,"#54ADF1"],[0.97324414715719054,"#54ADF2"],[0.97658862876254171,"#54AEF2"],[0.97993311036789288,"#55AEF3"],[0.98327759197324405,"#55AFF4"],[0.98662207357859522,"#55AFF4"],[0.98996655518394649,"#55B0F5"],[0.99331103678929766,"#56B0F6"],[0.99665551839464883,"#56B1F6"],[1,"#56B1F7"]],"colorbar":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"thickness":23.040000000000003,"title":"colData(se.subset)$TIN_median","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"tickmode":"array","ticktext":["50","60","70","80"],"tickvals":[0.11292248062015496,0.34470542635658902,0.57648837209302317,0.80827131782945727],"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"ticklen":2,"len":0.5,"yanchor":"top","y":1}},"xaxis":"x","yaxis":"y","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":37.260273972602747,"l":48.949771689497723},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.37622272276485136,0.39151493073170307],"tickmode":"array","ticktext":["-0.2","0.0","0.2"],"tickvals":[-0.20000000000000001,0,0.20000000000000007],"categoryorder":"array","categoryarray":["-0.2","0.0","0.2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":{"text":"PC1: 23% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.59358348614707002,0.53505418657770532],"tickmode":"array","ticktext":["-0.3","0.0","0.3"],"tickvals":[-0.30000000000000004,0,0.30000000000000004],"categoryorder":"array","categoryarray":["-0.3","0.0","0.3"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"PC2: 16% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"y":0.5,"yanchor":"top","title":{"text":"Adequacy_group<br />colData(se.subset)$TIN_median","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"14fc23038a03f":{"x":{},"y":{},"colour":{},"label":{},"shape":{},"type":"scatter"}},"cur_data":"14fc23038a03f","visdat":{"14fc23038a03f":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```


The PCA plot shows the first two principal components of the log cpm transformed expression data 
TIN values illustrates for the variation that is observed and shape indicates where samples originate. 

### PCA plot colored by RIN 


```r
plot3 <- ggplot(as.data.frame(pca$rotation[,c(1,2)]),
                aes(PC1,PC2, color=as.numeric(colData(se.subset)$RIN),label = sample_id)) +
  (if (!is.null(shape_var)) geom_point(size=3, aes(shape=colData(se.subset)[[shape_var]])) else geom_point(size=3)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  labs(shape = shape_var)
safe_ggplotly(plot3)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-9b3b48d2f01e6865b9ec" style="width:1056px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-9b3b48d2f01e6865b9ec">{"x":{"data":[{"x":[-0.34132555669682618,0.19930779475461491,-0.31978116431201509,0.26440980228685385,-0.21468381307301776,0.26409980068254812],"y":[0.002143706201169712,-0.54228177375048936,-0.013105187639462749,0.14215137559030044,0.4837524741811246,-0.16606817341186514],"text":["PC1: -0.3413256<br />PC2:  0.002143706<br />as.numeric(colData(se.subset)$RIN): 5.8<br />sample_id: Data2Cohort_P01_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1:  0.1993078<br />PC2: -0.542281774<br />as.numeric(colData(se.subset)$RIN): 8.0<br />sample_id: Data2Cohort_P02_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1: -0.3197812<br />PC2: -0.013105188<br />as.numeric(colData(se.subset)$RIN): 8.6<br />sample_id: Data2Cohort_P05_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1:  0.2644098<br />PC2:  0.142151376<br />as.numeric(colData(se.subset)$RIN): 8.7<br />sample_id: Data2Cohort_P06_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1: -0.2146838<br />PC2:  0.483752474<br />as.numeric(colData(se.subset)$RIN): 8.9<br />sample_id: Data2Cohort_P09_S1<br />colData(se.subset)[[shape_var]]: Adequate","PC1:  0.2640998<br />PC2: -0.166068173<br />as.numeric(colData(se.subset)$RIN): 7.5<br />sample_id: Data2Cohort_P10_S1<br />colData(se.subset)[[shape_var]]: Adequate"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(19,43,67,1)","rgba(65,135,191,1)","rgba(79,163,228,1)","rgba(81,168,234,1)","rgba(86,177,247,1)","rgba(54,113,161,1)"],"opacity":1,"size":11.338582677165356,"symbol":"circle","line":{"width":1.8897637795275593,"color":["rgba(19,43,67,1)","rgba(65,135,191,1)","rgba(79,163,228,1)","rgba(81,168,234,1)","rgba(86,177,247,1)","rgba(54,113,161,1)"]}},"hoveron":"points","name":"Adequate","legendgroup":"Adequate","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.33147941528946495,-0.19479135720787935,-0.29838270612673901],"y":[-0.16292302179947601,0.36308795005655825,-0.41814741950596496],"text":["PC1: -0.3314794<br />PC2: -0.162923022<br />as.numeric(colData(se.subset)$RIN): 5.9<br />sample_id: Data2Cohort_P03_S1<br />colData(se.subset)[[shape_var]]: Borderline","PC1: -0.1947914<br />PC2:  0.363087950<br />as.numeric(colData(se.subset)$RIN): 6.2<br />sample_id: Data2Cohort_P07_S1<br />colData(se.subset)[[shape_var]]: Borderline","PC1: -0.2983827<br />PC2: -0.418147420<br />as.numeric(colData(se.subset)$RIN): 8.1<br />sample_id: Data2Cohort_P11_S1<br />colData(se.subset)[[shape_var]]: Borderline"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(21,47,72,1)","rgba(27,58,88,1)","rgba(67,140,197,1)"],"opacity":1,"size":11.338582677165356,"symbol":"triangle-up","line":{"width":1.8897637795275593,"color":["rgba(21,47,72,1)","rgba(27,58,88,1)","rgba(67,140,197,1)"]}},"hoveron":"points","name":"Borderline","legendgroup":"Borderline","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33961684010506166,0.2763920102131891,0.35661776466367789],"y":[-0.043233468011760666,0.063338233143606607,0.29128530494625321],"text":["PC1:  0.3396168<br />PC2: -0.043233468<br />as.numeric(colData(se.subset)$RIN): 7.5<br />sample_id: Data2Cohort_P04_S1<br />colData(se.subset)[[shape_var]]: Intermediate","PC1:  0.2763920<br />PC2:  0.063338233<br />as.numeric(colData(se.subset)$RIN): 6.1<br />sample_id: Data2Cohort_P08_S1<br />colData(se.subset)[[shape_var]]: Intermediate","PC1:  0.3566178<br />PC2:  0.291285305<br />as.numeric(colData(se.subset)$RIN): 7.3<br />sample_id: Data2Cohort_P12_S1<br />colData(se.subset)[[shape_var]]: Intermediate"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(54,113,161,1)","rgba(25,54,82,1)","rgba(50,104,149,1)"],"opacity":1,"size":11.338582677165356,"symbol":"square","line":{"width":1.8897637795275593,"color":["rgba(54,113,161,1)","rgba(25,54,82,1)","rgba(50,104,149,1)"]}},"hoveron":"points","name":"Intermediate","legendgroup":"Intermediate","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.20000000000000001],"y":[-0.30000000000000004],"name":"0e348c48abe37d6f90493d1d5a22826d","type":"scatter","mode":"markers","opacity":0,"hoverinfo":"skip","showlegend":false,"marker":{"color":[0,1],"colorscale":[[0,"#132B43"],[0.003344481605351046,"#132B44"],[0.0066889632107023783,"#132C44"],[0.010033444816053425,"#142C45"],[0.013377926421404757,"#142D45"],[0.016722408026755804,"#142D46"],[0.020066889632107135,"#142D46"],[0.023411371237458182,"#142E47"],[0.026755852842809229,"#152E47"],[0.03010033444816056,"#152F48"],[0.033444816053511607,"#152F48"],[0.036789297658862942,"#152F49"],[0.040133779264213985,"#153049"],[0.04347826086956532,"#16304A"],[0.046822742474916364,"#16304A"],[0.050167224080267699,"#16314B"],[0.053511705685618742,"#16314B"],[0.056856187290969785,"#16324C"],[0.06020066889632112,"#17324D"],[0.063545150501672171,"#17324D"],[0.066889632107023492,"#17334E"],[0.070234113712374549,"#17334E"],[0.073578595317725884,"#17344F"],[0.076923076923076927,"#18344F"],[0.080267558528427971,"#183450"],[0.083612040133779306,"#183550"],[0.086956521739130349,"#183551"],[0.090301003344481684,"#183651"],[0.093645484949832727,"#193652"],[0.096989966555184062,"#193652"],[0.10033444816053511,"#193753"],[0.10367892976588644,"#193754"],[0.10702341137123748,"#193854"],[0.11036789297658853,"#1A3855"],[0.11371237458193986,"#1A3955"],[0.11705685618729091,"#1A3956"],[0.12040133779264224,"#1A3956"],[0.12374581939799328,"#1A3A57"],[0.12709030100334462,"#1B3A57"],[0.13043478260869568,"#1B3B58"],[0.13377926421404671,"#1B3B59"],[0.13712374581939804,"#1B3B59"],[0.1404682274247491,"#1C3C5A"],[0.14381270903010043,"#1C3C5A"],[0.14715719063545146,"#1C3D5B"],[0.1505016722408028,"#1C3D5B"],[0.15384615384615385,"#1C3D5C"],[0.15719063545150488,"#1D3E5C"],[0.16053511705685622,"#1D3E5D"],[0.16387959866220728,"#1D3F5D"],[0.16722408026755861,"#1D3F5E"],[0.17056856187290964,"#1D3F5F"],[0.17391304347826098,"#1E405F"],[0.17725752508361203,"#1E4060"],[0.18060200668896306,"#1E4160"],[0.1839464882943144,"#1E4161"],[0.18729096989966545,"#1E4261"],[0.19063545150501679,"#1F4262"],[0.19397993311036782,"#1F4263"],[0.19732441471571915,"#1F4363"],[0.20066889632107021,"#1F4364"],[0.20401337792642155,"#1F4464"],[0.20735785953177258,"#204465"],[0.21070234113712391,"#204465"],[0.21404682274247497,"#204566"],[0.21739130434782603,"#204566"],[0.22073578595317733,"#214667"],[0.22408026755852839,"#214668"],[0.22742474916387972,"#214768"],[0.23076923076923078,"#214769"],[0.23411371237458209,"#214769"],[0.23745819397993315,"#22486A"],[0.2408026755852842,"#22486A"],[0.24414715719063554,"#22496B"],[0.24749163879598657,"#22496C"],[0.2508361204013379,"#224A6C"],[0.25418060200668896,"#234A6D"],[0.2575250836120403,"#234A6D"],[0.26086956521739135,"#234B6E"],[0.26421404682274235,"#234B6E"],[0.26755852842809369,"#244C6F"],[0.27090301003344475,"#244C70"],[0.27424749163879608,"#244C70"],[0.27759197324414714,"#244D71"],[0.28093645484949847,"#244D71"],[0.28428093645484953,"#254E72"],[0.28762541806020053,"#254E72"],[0.29096989966555187,"#254F73"],[0.29431438127090293,"#254F74"],[0.29765886287625426,"#254F74"],[0.30100334448160532,"#265075"],[0.30434782608695665,"#265075"],[0.30769230769230771,"#265176"],[0.31103678929765871,"#265176"],[0.31438127090301005,"#275277"],[0.3177257525083611,"#275278"],[0.32107023411371244,"#275278"],[0.3244147157190635,"#275379"],[0.32775919732441483,"#275379"],[0.33110367892976589,"#28547A"],[0.33444816053511695,"#28547B"],[0.33779264214046828,"#28557B"],[0.34113712374581956,"#28557C"],[0.34448160535117062,"#28567C"],[0.34782608695652167,"#29567D"],[0.35117056856187301,"#29567D"],[0.35451505016722407,"#29577E"],[0.35785953177257512,"#29577F"],[0.36120401337792646,"#2A587F"],[0.36454849498327779,"#2A5880"],[0.3678929765886288,"#2A5980"],[0.37123745819397985,"#2A5981"],[0.37458193979933119,"#2A5982"],[0.37792642140468224,"#2B5A82"],[0.38127090301003358,"#2B5A83"],[0.38461538461538464,"#2B5B83"],[0.38795986622073597,"#2B5B84"],[0.39130434782608697,"#2C5C85"],[0.39464882943143803,"#2C5C85"],[0.39799331103678937,"#2C5D86"],[0.40133779264214042,"#2C5D86"],[0.40468227424749176,"#2C5D87"],[0.40802675585284282,"#2D5E87"],[0.41137123745819415,"#2D5E88"],[0.41471571906354515,"#2D5F89"],[0.41806020066889621,"#2D5F89"],[0.42140468227424754,"#2E608A"],[0.4247491638795986,"#2E608A"],[0.42809364548494994,"#2E618B"],[0.43143812709030099,"#2E618C"],[0.43478260869565233,"#2E618C"],[0.43812709030100333,"#2F628D"],[0.44147157190635439,"#2F628D"],[0.44481605351170572,"#2F638E"],[0.44816053511705678,"#2F638F"],[0.45150501672240811,"#30648F"],[0.45484949832775917,"#306490"],[0.45819397993311051,"#306590"],[0.46153846153846156,"#306591"],[0.46488294314381257,"#306592"],[0.4682274247491639,"#316692"],[0.47157190635451524,"#316693"],[0.47491638795986629,"#316793"],[0.47826086956521735,"#316794"],[0.48160535117056869,"#326895"],[0.48494983277591974,"#326895"],[0.48829431438127074,"#326996"],[0.49163879598662208,"#326996"],[0.49498327759197341,"#326997"],[0.49832775919732447,"#336A98"],[0.50167224080267547,"#336A98"],[0.50501672240802686,"#336B99"],[0.50836120401337792,"#336B99"],[0.5117056856187292,"#346C9A"],[0.51505016722408026,"#346C9B"],[0.51839464882943165,"#346D9B"],[0.52173913043478271,"#346D9C"],[0.52508361204013365,"#346E9D"],[0.52842809364548504,"#356E9D"],[0.5317725752508361,"#356E9E"],[0.53511705685618738,"#356F9E"],[0.53846153846153844,"#356F9F"],[0.54180602006688983,"#3670A0"],[0.54515050167224088,"#3670A0"],[0.54849498327759194,"#3671A1"],[0.55183946488294322,"#3671A1"],[0.55518394648829428,"#3772A2"],[0.55852842809364556,"#3772A3"],[0.56187290969899661,"#3773A3"],[0.565217391304348,"#3773A4"],[0.56856187290969906,"#3773A4"],[0.57190635451505012,"#3874A5"],[0.5752508361204014,"#3874A6"],[0.57859531772575246,"#3875A6"],[0.58193979933110374,"#3875A7"],[0.58528428093645479,"#3976A8"],[0.58862876254180618,"#3976A8"],[0.59197324414715724,"#3977A9"],[0.5953177257525083,"#3977A9"],[0.59866220735785958,"#3978AA"],[0.60200668896321097,"#3A78AB"],[0.60535117056856191,"#3A79AB"],[0.60869565217391297,"#3A79AC"],[0.61204013377926436,"#3A79AC"],[0.61538461538461542,"#3B7AAD"],[0.61872909698996648,"#3B7AAE"],[0.62207357859531776,"#3B7BAE"],[0.62541806020066915,"#3B7BAF"],[0.62876254180602009,"#3C7CB0"],[0.63210702341137115,"#3C7CB0"],[0.63545150501672254,"#3C7DB1"],[0.6387959866220736,"#3C7DB1"],[0.64214046822742465,"#3C7EB2"],[0.64548494983277593,"#3D7EB3"],[0.64882943143812732,"#3D7FB3"],[0.65217391304347827,"#3D7FB4"],[0.65551839464882933,"#3D7FB5"],[0.65886287625418072,"#3E80B5"],[0.66220735785953178,"#3E80B6"],[0.66555183946488283,"#3E81B6"],[0.66889632107023411,"#3E81B7"],[0.6722408026755855,"#3F82B8"],[0.67558528428093656,"#3F82B8"],[0.67892976588628751,"#3F83B9"],[0.6822742474916389,"#3F83BA"],[0.68561872909699018,"#4084BA"],[0.68896321070234101,"#4084BB"],[0.69230769230769229,"#4085BB"],[0.69565217391304368,"#4085BC"],[0.69899665551839474,"#4086BD"],[0.70234113712374568,"#4186BD"],[0.70568561872909707,"#4186BE"],[0.70903010033444835,"#4187BF"],[0.71237458193979974,"#4187BF"],[0.71571906354515047,"#4288C0"],[0.71906354515050186,"#4288C1"],[0.72240802675585258,"#4289C1"],[0.72575250836120386,"#4289C2"],[0.72909698996655525,"#438AC2"],[0.73244147157190653,"#438AC3"],[0.73578595317725792,"#438BC4"],[0.73913043478260865,"#438BC4"],[0.74247491638796004,"#438CC5"],[0.74581939799331076,"#448CC6"],[0.74916387959866204,"#448DC6"],[0.75250836120401343,"#448DC7"],[0.75585284280936471,"#448EC8"],[0.7591973244147161,"#458EC8"],[0.76254180602006683,"#458FC9"],[0.76588628762541822,"#458FC9"],[0.76923076923076894,"#458FCA"],[0.77257525083612033,"#4690CB"],[0.77591973244147161,"#4690CB"],[0.77926421404682289,"#4691CC"],[0.78260869565217428,"#4691CD"],[0.785953177257525,"#4792CD"],[0.78929765886287639,"#4792CE"],[0.79264214046822712,"#4793CF"],[0.79598662207357851,"#4793CF"],[0.79933110367892979,"#4894D0"],[0.80267558528428107,"#4894D0"],[0.80602006688963246,"#4895D1"],[0.80936454849498318,"#4895D2"],[0.81270903010033457,"#4896D2"],[0.81605351170568585,"#4996D3"],[0.81939799331103669,"#4997D4"],[0.82274247491638797,"#4997D4"],[0.82608695652173936,"#4998D5"],[0.82943143812709064,"#4A98D6"],[0.83277591973244136,"#4A99D6"],[0.83612040133779275,"#4A99D7"],[0.83946488294314403,"#4A9AD8"],[0.84280936454849487,"#4B9AD8"],[0.84615384615384615,"#4B9BD9"],[0.84949832775919754,"#4B9BDA"],[0.85284280936454882,"#4B9BDA"],[0.85618729096989954,"#4C9CDB"],[0.85953177257525093,"#4C9CDB"],[0.86287625418060221,"#4C9DDC"],[0.86622073578595304,"#4C9DDD"],[0.86956521739130432,"#4D9EDD"],[0.87290969899665571,"#4D9EDE"],[0.87625418060200699,"#4D9FDF"],[0.87959866220735772,"#4D9FDF"],[0.88294314381270911,"#4DA0E0"],[0.88628762541806039,"#4EA0E1"],[0.88963210702341122,"#4EA1E1"],[0.8929765886287625,"#4EA1E2"],[0.89632107023411389,"#4EA2E3"],[0.89966555183946517,"#4FA2E3"],[0.9030100334448159,"#4FA3E4"],[0.90635451505016729,"#4FA3E5"],[0.90969899665551857,"#4FA4E5"],[0.9130434782608694,"#50A4E6"],[0.91638795986622068,"#50A5E7"],[0.91973244147157207,"#50A5E7"],[0.92307692307692335,"#50A6E8"],[0.92642140468227407,"#51A6E8"],[0.92976588628762546,"#51A7E9"],[0.93311036789297674,"#51A7EA"],[0.93645484949832758,"#51A8EA"],[0.93979933110367886,"#52A8EB"],[0.94314381270903025,"#52A9EC"],[0.94648829431438153,"#52A9EC"],[0.94983277591973225,"#52AAED"],[0.95317725752508364,"#53AAEE"],[0.95652173913043492,"#53ABEE"],[0.95986622073578576,"#53ABEF"],[0.96321070234113704,"#53ACF0"],[0.96655518394648843,"#54ACF0"],[0.96989966555183971,"#54ADF1"],[0.97324414715719043,"#54ADF2"],[0.97658862876254182,"#54AEF2"],[0.9799331103678931,"#55AEF3"],[0.98327759197324394,"#55AFF4"],[0.98662207357859522,"#55AFF4"],[0.98996655518394661,"#55B0F5"],[0.99331103678929789,"#56B0F6"],[0.99665551839464928,"#56B1F6"],[1,"#56B1F7"]],"colorbar":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"thickness":23.040000000000003,"title":"as.numeric(colData(se.subset)$RIN)","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"tickmode":"array","ticktext":["6","7","8"],"tickvals":[0.065967741935483909,0.38747311827956982,0.70897849462365581],"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"ticklen":2,"len":0.5,"yanchor":"top","y":1}},"xaxis":"x","yaxis":"y","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":37.260273972602747,"l":48.949771689497723},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.37622272276485136,0.39151493073170307],"tickmode":"array","ticktext":["-0.2","0.0","0.2"],"tickvals":[-0.20000000000000001,0,0.20000000000000007],"categoryorder":"array","categoryarray":["-0.2","0.0","0.2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":{"text":"PC1: 23% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.59358348614707002,0.53505418657770532],"tickmode":"array","ticktext":["-0.3","0.0","0.3"],"tickvals":[-0.30000000000000004,0,0.30000000000000004],"categoryorder":"array","categoryarray":["-0.3","0.0","0.3"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"PC2: 16% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"y":0.5,"yanchor":"top","title":{"text":"as.numeric(colData(se.subset)$RIN)<br />Adequacy_group","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"14fc239761b9a":{"x":{},"y":{},"colour":{},"label":{},"shape":{},"type":"scatter"}},"cur_data":"14fc239761b9a","visdat":{"14fc239761b9a":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```


The PCA plot shows the first two principal components of the log cpm transformed expression data 
RIN values illustrates for the variation that is observed and shape indicates where samples originate. 

### PCA plot colored by outliers 


```r
plot4 <- ggplot(as.data.frame(pca$rotation[,c(1,2)]),
                aes(PC1,PC2, label = sample_id)) +
  (if (!is.null(shape_var)) aes(shape = colData(se.subset)[[shape_var]]) else NULL) +
  geom_point(size=4, aes(color = out.color)) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 14)
  ) +
  ggplot2::xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ggplot2::ylab(paste0("PC2 (", percentVar[2], "%)")) +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 15), legend.text = ggplot2::element_text(size = 13)) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2))) +
  labs(color='Outliers', shape = shape_var)+
  ggplot2::scale_color_manual(name = "Outliers",
                     values = c("firebrick", "blue"),
                     labels = c("no", "yes"))


safe_ggplotly(plot4)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-7ec0d789fca32ee51100" style="width:960px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-7ec0d789fca32ee51100">{"x":{"data":[{"x":[-0.34132555669682618,0.19930779475461491,-0.31978116431201509,0.26440980228685385,-0.21468381307301776,0.26409980068254812],"y":[0.002143706201169712,-0.54228177375048936,-0.013105187639462749,0.14215137559030044,0.4837524741811246,-0.16606817341186514],"text":["colData(se.subset)[[shape_var]]: Adequate<br />PC1: -0.3413256<br />PC2:  0.002143706<br />sample_id: Data2Cohort_P01_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Adequate<br />PC1:  0.1993078<br />PC2: -0.542281774<br />sample_id: Data2Cohort_P02_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Adequate<br />PC1: -0.3197812<br />PC2: -0.013105188<br />sample_id: Data2Cohort_P05_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Adequate<br />PC1:  0.2644098<br />PC2:  0.142151376<br />sample_id: Data2Cohort_P06_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Adequate<br />PC1: -0.2146838<br />PC2:  0.483752474<br />sample_id: Data2Cohort_P09_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Adequate<br />PC1:  0.2640998<br />PC2: -0.166068173<br />sample_id: Data2Cohort_P10_S1<br />out.color: no"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(178,34,34,1)","opacity":1,"size":15.118110236220474,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(178,34,34,1)"}},"hoveron":"points","name":"(no,Adequate)","legendgroup":"(no,Adequate)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.33147941528946495,-0.19479135720787935,-0.29838270612673901],"y":[-0.16292302179947601,0.36308795005655825,-0.41814741950596496],"text":["colData(se.subset)[[shape_var]]: Borderline<br />PC1: -0.3314794<br />PC2: -0.162923022<br />sample_id: Data2Cohort_P03_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Borderline<br />PC1: -0.1947914<br />PC2:  0.363087950<br />sample_id: Data2Cohort_P07_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Borderline<br />PC1: -0.2983827<br />PC2: -0.418147420<br />sample_id: Data2Cohort_P11_S1<br />out.color: no"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(178,34,34,1)","opacity":1,"size":15.118110236220474,"symbol":"triangle-up","line":{"width":1.8897637795275593,"color":"rgba(178,34,34,1)"}},"hoveron":"points","name":"(no,Borderline)","legendgroup":"(no,Borderline)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.33961684010506166,0.2763920102131891,0.35661776466367789],"y":[-0.043233468011760666,0.063338233143606607,0.29128530494625321],"text":["colData(se.subset)[[shape_var]]: Intermediate<br />PC1:  0.3396168<br />PC2: -0.043233468<br />sample_id: Data2Cohort_P04_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Intermediate<br />PC1:  0.2763920<br />PC2:  0.063338233<br />sample_id: Data2Cohort_P08_S1<br />out.color: no","colData(se.subset)[[shape_var]]: Intermediate<br />PC1:  0.3566178<br />PC2:  0.291285305<br />sample_id: Data2Cohort_P12_S1<br />out.color: no"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(178,34,34,1)","opacity":1,"size":15.118110236220474,"symbol":"square","line":{"width":1.8897637795275593,"color":"rgba(178,34,34,1)"}},"hoveron":"points","name":"(no,Intermediate)","legendgroup":"(no,Intermediate)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":48.152760481527615,"l":66.749688667496883},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.37622272276485136,0.39151493073170307],"tickmode":"array","ticktext":["-0.2","0.0","0.2"],"tickvals":[-0.20000000000000001,0,0.20000000000000007],"categoryorder":"array","categoryarray":["-0.2","0.0","0.2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":18.596928185969286},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"PC1 (23%)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.59358348614707002,0.53505418657770532],"tickmode":"array","ticktext":["-0.3","0.0","0.3"],"tickvals":[-0.30000000000000004,0,0.30000000000000004],"categoryorder":"array","categoryarray":["-0.3","0.0","0.3"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":18.596928185969279},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"PC2 (16%)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"rgba(255,255,255,1)","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":17.268576172685766},"title":{"text":"Adequacy_group<br />Outliers","font":{"color":"rgba(0,0,0,1)","family":"","size":19.925280199252807}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"14fc278e88a3e":{"shape":{},"x":{},"y":{},"label":{},"colour":{},"type":"scatter"}},"cur_data":"14fc278e88a3e","visdat":{"14fc278e88a3e":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

The PCA plot shows the first two principal components of the log cpm transformed expression data 
Color illustrates if samples are outliers based on -+2 standard deviation away from the mean and shape indicates where samples originate. 

## Batch information {.tabset .tabset-fade}
### PCA batch  
color PCA plot by batch and see if you observe visual separation of the samples 

```r
batch.sub <- as.factor(colData(se.subset)$batch_id)
nb <- length(unique(batch.sub))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb)
myCol <- setNames(mycolors,unique(batch.sub))

pca.before <- prcomp(z.data,scale=F, center=F)

perc.var <- (pca.before$sdev)^2/sum(pca.before$sdev^2)
percentVar <- round(100 * perc.var)
sampleID <- rownames(pca.before$rotation)
# The tissue_type shape here was always single-level (se.subset is already
# tissue-filtered, so it's purely decorative even when present) -- drop it entirely
# when tissue_type isn't available at all.
batch.plot <- ggplot(as.data.frame(pca.before$rotation[,c(1,2)]),
                     aes(PC1,PC2, color=batch.sub, label=sampleID)) +
  (if (has_tissue_type) aes(shape = colData(se.subset)$tissue_type) else NULL) +
  geom_point(size=3) +
  labs(#title = "Before batch effect correction",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance")) +
   (if (has_tissue_type) scale_shape_manual(name = "Sample Type", values = 16, labels = tissue) else NULL) +
  scale_color_manual(name = if (is.null(tissue)) tissueType else tissue,
                     #values = COLORS
                     values = myCol)+
  theme_classic() +
  theme(text = element_text(size = 14))


safe_ggplotly(batch.plot)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-005ea33c2266a56e9e6d" style="width:960px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-005ea33c2266a56e9e6d">{"x":{"data":[{"x":[-0.34132555669682618,-0.33147941528946495,-0.31978116431201509,-0.19479135720787935,-0.21468381307301776,-0.29838270612673901],"y":[0.002143706201169712,-0.16292302179947601,-0.013105187639462749,0.36308795005655825,0.4837524741811246,-0.41814741950596496],"text":["PC1: -0.3413256<br />PC2:  0.002143706<br />batch.sub: Batch1<br />sampleID: Data2Cohort_P01_S1","PC1: -0.3314794<br />PC2: -0.162923022<br />batch.sub: Batch1<br />sampleID: Data2Cohort_P03_S1","PC1: -0.3197812<br />PC2: -0.013105188<br />batch.sub: Batch1<br />sampleID: Data2Cohort_P05_S1","PC1: -0.1947914<br />PC2:  0.363087950<br />batch.sub: Batch1<br />sampleID: Data2Cohort_P07_S1","PC1: -0.2146838<br />PC2:  0.483752474<br />batch.sub: Batch1<br />sampleID: Data2Cohort_P09_S1","PC1: -0.2983827<br />PC2: -0.418147420<br />batch.sub: Batch1<br />sampleID: Data2Cohort_P11_S1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(166,206,227,1)","opacity":1,"size":11.338582677165356,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(166,206,227,1)"}},"hoveron":"points","name":"Batch1","legendgroup":"Batch1","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.19930779475461491,0.33961684010506166,0.26440980228685385,0.2763920102131891,0.26409980068254812,0.35661776466367789],"y":[-0.54228177375048936,-0.043233468011760666,0.14215137559030044,0.063338233143606607,-0.16606817341186514,0.29128530494625321],"text":["PC1:  0.1993078<br />PC2: -0.542281774<br />batch.sub: Batch2<br />sampleID: Data2Cohort_P02_S1","PC1:  0.3396168<br />PC2: -0.043233468<br />batch.sub: Batch2<br />sampleID: Data2Cohort_P04_S1","PC1:  0.2644098<br />PC2:  0.142151376<br />batch.sub: Batch2<br />sampleID: Data2Cohort_P06_S1","PC1:  0.2763920<br />PC2:  0.063338233<br />batch.sub: Batch2<br />sampleID: Data2Cohort_P08_S1","PC1:  0.2640998<br />PC2: -0.166068173<br />batch.sub: Batch2<br />sampleID: Data2Cohort_P10_S1","PC1:  0.3566178<br />PC2:  0.291285305<br />batch.sub: Batch2<br />sampleID: Data2Cohort_P12_S1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(177,89,40,1)","opacity":1,"size":11.338582677165356,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(177,89,40,1)"}},"hoveron":"points","name":"Batch2","legendgroup":"Batch2","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":44.433374844333756,"l":59.310917393109179},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.37622272276485136,0.39151493073170307],"tickmode":"array","ticktext":["-0.2","0.0","0.2"],"tickvals":[-0.20000000000000001,0,0.20000000000000007],"categoryorder":"array","categoryarray":["-0.2","0.0","0.2"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":14.87754254877543},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176002,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"PC1: 23% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.59358348614707002,0.53505418657770532],"tickmode":"array","ticktext":["-0.3","0.0","0.3"],"tickvals":[-0.30000000000000004,0,0.30000000000000004],"categoryorder":"array","categoryarray":["-0.3","0.0","0.3"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":14.877542548775427},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176002,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"PC2: 16% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"shapes":[],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":14.87754254877543},"title":{"text":"Data2-style Synthetic Example QC","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"14fc23debf921":{"x":{},"y":{},"colour":{},"label":{},"type":"scatter"}},"cur_data":"14fc23debf921","visdat":{"14fc23debf921":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

## QC Tally and Number of Outliers Binned {.tabset .tabset-fade}

### QC histogram outliers 


```r
# ------------------------------------
# 1. Remove QC Metrics with Constant Values
# ------------------------------------

# Identify columns with the same values across all samples
rm.ix.qc <- c()
for (i in 1:ncol(qcMetrics.subset)) {
  if (length(unique(qcMetrics.subset[,i]))== 1) {
    rm.ix.qc <- c(rm.ix.qc, i)
  }
}
# Filter dataset to remove constant QC metrics
# NOTE: `-rm.ix.qc` errors ("invalid argument to unary operator") when rm.ix.qc is
# NULL, i.e. when no QC metric happens to be constant across this tissue subset
# -- guard so that case (a real possibility, not just an edge case) doesn't crash.
qc.metrics.unique.col <- if (length(rm.ix.qc) > 0) qcMetrics.subset[-rm.ix.qc] else qcMetrics.subset

# ------------------------------------
# 2. Identify Outliers in QC Metrics
# ------------------------------------

# Initialize matrices for tallying outliers
qc.tally <- matrix(ncol = nrow(qc.metrics.unique.col), data = 0)
colnames(qc.tally) <- rownames(qc.metrics.unique.col)

qc.tally2 <- matrix(
  ncol = nrow(qc.metrics.unique.col),
  nrow = ncol(qc.metrics.unique.col),
  data = 0)
rownames(qc.tally2) <- colnames(qc.metrics.unique.col)
colnames(qc.tally2) <- rownames(qc.metrics.unique.col)

isNum <- 0
# Loop through each QC metric column
for (i in 1:ncol(qc.metrics.unique.col)){
  qc <- qc.metrics.unique.col[,i]
  if (class(qc) == "numeric" | class(qc) == "integer") {
    isNum <- isNum + 1
    mean.qc <- mean(qc)
    sd.qc <- sd(qc)
    # Identify samples ±2 SD from mean correlation
    outlier_samples <- which(qc > mean.qc + 2 * sd.qc | qc < mean.qc - 2 * sd.qc)
    # Increment tally matrices
    qc.tally[,outlier_samples] <- qc.tally[,outlier_samples] + 1
    qc.tally2[i, outlier_samples] <- qc.tally2[i, outlier_samples] + 1
  }
}
# ------------------------------------
# 3. Prepare Data for Histogram Plot
# ------------------------------------

# Convert tally matrix into a data frame
qc.tally.gg <- as.data.frame(t(qc.tally))
# Create a column to indicate outlier binning
qc.tally.gg$Number_Outliers_Binned <- qc.tally.gg$V1 > 6
qc.tally.gg$FillColor <- sapply(qc.tally.gg$Number_Outliers_Binned, function(x) {
  if(x == TRUE) {
    return("Above 6 Outliers")
  } else{
    return("Below 6 Outliers")
  }
})
# ------------------------------------
# 4. Plot Histogram of Outlier Tally
# ------------------------------------
histQC <- ggplot(data = qc.tally.gg, aes(x = V1, fill = FillColor, color = FillColor)) +
  geom_histogram() +
  labs(x = "Tally of outliers") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 15),
    axis.title = ggplot2::element_text(size = 15),
    legend.title = element_blank()
  )

histQC
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/histo_tally_outliers-1.png" style="display: block; margin: auto;" />

```r
#ggplotly(histQC)
```

The plot indicates tally of QC metrics for each sample. It shows a separation where number of outliers are binned as below 6 (green) or above 6 (red).  



```r
qc.tally.gg$Number_Outliers_Binned <- as.character(qc.tally.gg$Number_Outliers_Binned)
qc.tally.gg$Number_Outliers_Binned <- gsub("FALSE", "x < 6", qc.tally.gg$Number_Outliers_Binned)
qc.tally.gg$Number_Outliers_Binned <- gsub("TRUE", "x => 6", qc.tally.gg$Number_Outliers_Binned)

highAmountOutlier <- colnames(qc.tally)[which(qc.tally > 6)]

qualAnnot <- sapply(qc.tally, function(x){
  if ( x < 6){
    return("x =< 6")
  }else if (x < 40){
    return("6 < x < 20")
  }else {
    return("x >= 20")
  }
})

Number_Outliers_Binned <- as.factor(qualAnnot)
Number_Outliers_Binned <- factor(Number_Outliers_Binned, 
                                 levels = c("x =< 6", 
                                            "6 < x < 20", 
                                            "x >= 20"))
```


A total of 26 QC metrics were taken into consideration. The following samples had a high number of QC metrics outliers (more than 6): . These are colored as red on the histogram.


### Correlaton of uniquely mapped percent  


```r
qc.df2 <- as_tibble(qcMetrics.subset@listData)
qc.df2$Number_Outliers_Binned <- Number_Outliers_Binned
ggscatter(qc.df2, x = "STAR_uniquely_mapped_percent", y = "rnaseqc_Unique.Rate.of.Mapped", color = "Number_Outliers_Binned",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "STAR uniquely mapped percent", ylab = "RNASeQC unique percent")
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/uniquely_mapped_percent-1.png" style="display: block; margin: auto;" />

Correlation between RSeQC unique mapped percent and STAR uniquely mapped pecent  colored by number of outliers binned.


### Correlaton of RIN vs TIN 


```r
ggscatter(qc.df2, x = "RIN", y = "TIN_median", color = "Number_Outliers_Binned",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "RNA integrity number (RIN)", ylab = "Transcript integrity number (TIN")
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/RIN_vs_TIN_median-1.png" style="display: block; margin: auto;" />

Correlation between Transcript integrity number (TIN) and RNA integrity number (RIN) colored by number of outliers binned.
 

### Correlation uniquelly mapped percent vs genes detected 


```r
qc.df2$Genes_detected <- colSums(assay(se.subset) > 0)
ggscatter(qc.df2, x = "STAR_uniquely_mapped_percent", y = "Genes_detected", color = "Number_Outliers_Binned",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "STAR uniquely mapped percent", ylab = "Genes detected")
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/uniquely_mapped_percent_vs_gene_detected-1.png" style="display: block; margin: auto;" />

Correlation between Genes detected and STAR uniquely mapped percent colored by number of outliers binned.
Number of QC oultier binned facilities the identification of low quality samples possible duplicates as seen with two samples colored blue with high genes detected and low uniquely mapped percent.  

### PCA of QC metrics  


```r
#z-score scale and removing "SAMPLE_ID" column ( non-numeric) from qc.metrics 
qcMetrics.z <- t(scale(qcMetrics.subset))
#qcMetrics.z <- t(scale(qc.metrics))
# Eliminate any rows that have NA after normalization
qcMetrics.z <- qcMetrics.z[which(rowSums(qcMetrics.z) != "NaN"),]
qcMetrics.z.trim <- qcMetrics.z

qcMetrics.z.trim[qcMetrics.z < -2] <- -2
qcMetrics.z.trim[qcMetrics.z > 2] <- 2

Number_Outliers_Binned <- as.factor(qualAnnot)
Number_Outliers_Binned <- factor(Number_Outliers_Binned, 
                                 levels = c("x =< 6",
                                            "6 < x < 20",
                                            "x >= 20"))

z <- qcMetrics.z

#Run PCA
pc <- prcomp(z, center = F, scale=F)

pcabothtype <- as.data.frame(pc$rotation[,1:3])
pcabothtype$Number_Outliers_Binned <- Number_Outliers_Binned

spc <- summary(pc)

pc1Importance <- signif(spc$importance[2, 1] * 100, 3)
pc2Importance <- signif(spc$importance[2, 2] * 100, 3)
sampleID <- rownames(pcabothtype)
plotpca <- ggplot(pcabothtype, 
                  aes(x= PC1, y = PC2, color=as.numeric(colData(se.subset)$RIN),label = sampleID)) +
  geom_point(size=5, aes(shape=Number_Outliers_Binned))+
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 14),
    axis.title = ggplot2::element_text(size = 14)
  ) +
  ggplot2::xlab(paste0("PC1 (", pc1Importance, "%)")) +
  ggplot2::ylab(paste0("PC2 (", pc2Importance, "%)")) +
  ggplot2::theme(legend.title = ggplot2::element_text(size = 15), legend.text = ggplot2::element_text(size = 13)) #+
  #ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 2)))
#plotpca
safe_ggplotly(plotpca)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-b6e4b381ec93b71f4087" style="width:1152px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-b6e4b381ec93b71f4087">{"x":{"data":[{"x":[0.23002269327879263,0.33942068435323508,-0.13935469243904849,-0.075016190550066961,0.44109927877799776,0.30086549635366699,-0.22540927305826647,-0.45934681469319444,0.29906642267675365,-0.23003647583804307,-0.27024303225681107,-0.21106809660501466],"y":[0.35597781612936086,-0.019843377387354286,0.36356261232001846,-0.0056531774886370619,-0.12045803093960848,-0.48817228759425491,0.54005743206025236,-0.25703944800995443,0.17000403333185232,-0.23243563979166756,-0.16672947060474028,-0.13927046202526588],"text":["PC1:  0.23002269<br />PC2:  0.355977816<br />as.numeric(colData(se.subset)$RIN): 5.8<br />sampleID: Data2Cohort_P01_S1<br />Number_Outliers_Binned: x =< 6","PC1:  0.33942068<br />PC2: -0.019843377<br />as.numeric(colData(se.subset)$RIN): 8.0<br />sampleID: Data2Cohort_P02_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.13935469<br />PC2:  0.363562612<br />as.numeric(colData(se.subset)$RIN): 5.9<br />sampleID: Data2Cohort_P03_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.07501619<br />PC2: -0.005653177<br />as.numeric(colData(se.subset)$RIN): 7.5<br />sampleID: Data2Cohort_P04_S1<br />Number_Outliers_Binned: x =< 6","PC1:  0.44109928<br />PC2: -0.120458031<br />as.numeric(colData(se.subset)$RIN): 8.6<br />sampleID: Data2Cohort_P05_S1<br />Number_Outliers_Binned: x =< 6","PC1:  0.30086550<br />PC2: -0.488172288<br />as.numeric(colData(se.subset)$RIN): 8.7<br />sampleID: Data2Cohort_P06_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.22540927<br />PC2:  0.540057432<br />as.numeric(colData(se.subset)$RIN): 6.2<br />sampleID: Data2Cohort_P07_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.45934681<br />PC2: -0.257039448<br />as.numeric(colData(se.subset)$RIN): 6.1<br />sampleID: Data2Cohort_P08_S1<br />Number_Outliers_Binned: x =< 6","PC1:  0.29906642<br />PC2:  0.170004033<br />as.numeric(colData(se.subset)$RIN): 8.9<br />sampleID: Data2Cohort_P09_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.23003648<br />PC2: -0.232435640<br />as.numeric(colData(se.subset)$RIN): 7.5<br />sampleID: Data2Cohort_P10_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.27024303<br />PC2: -0.166729471<br />as.numeric(colData(se.subset)$RIN): 8.1<br />sampleID: Data2Cohort_P11_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.21106810<br />PC2: -0.139270462<br />as.numeric(colData(se.subset)$RIN): 7.3<br />sampleID: Data2Cohort_P12_S1<br />Number_Outliers_Binned: x =< 6"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(19,43,67,1)","rgba(65,135,191,1)","rgba(21,47,72,1)","rgba(54,113,161,1)","rgba(79,163,228,1)","rgba(81,168,234,1)","rgba(27,58,88,1)","rgba(25,54,82,1)","rgba(86,177,247,1)","rgba(54,113,161,1)","rgba(67,140,197,1)","rgba(50,104,149,1)"],"opacity":1,"size":18.897637795275593,"symbol":"circle","line":{"width":1.8897637795275593,"color":["rgba(19,43,67,1)","rgba(65,135,191,1)","rgba(21,47,72,1)","rgba(54,113,161,1)","rgba(79,163,228,1)","rgba(81,168,234,1)","rgba(27,58,88,1)","rgba(25,54,82,1)","rgba(86,177,247,1)","rgba(54,113,161,1)","rgba(67,140,197,1)","rgba(50,104,149,1)"]}},"hoveron":"points","name":"x =< 6","legendgroup":"x =< 6","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.5],"y":[-0.30000000000000004],"name":"0e348c48abe37d6f90493d1d5a22826d","type":"scatter","mode":"markers","opacity":0,"hoverinfo":"skip","showlegend":false,"marker":{"color":[0,1],"colorscale":[[0,"#132B43"],[0.003344481605351046,"#132B44"],[0.0066889632107023783,"#132C44"],[0.010033444816053425,"#142C45"],[0.013377926421404757,"#142D45"],[0.016722408026755804,"#142D46"],[0.020066889632107135,"#142D46"],[0.023411371237458182,"#142E47"],[0.026755852842809229,"#152E47"],[0.03010033444816056,"#152F48"],[0.033444816053511607,"#152F48"],[0.036789297658862942,"#152F49"],[0.040133779264213985,"#153049"],[0.04347826086956532,"#16304A"],[0.046822742474916364,"#16304A"],[0.050167224080267699,"#16314B"],[0.053511705685618742,"#16314B"],[0.056856187290969785,"#16324C"],[0.06020066889632112,"#17324D"],[0.063545150501672171,"#17324D"],[0.066889632107023492,"#17334E"],[0.070234113712374549,"#17334E"],[0.073578595317725884,"#17344F"],[0.076923076923076927,"#18344F"],[0.080267558528427971,"#183450"],[0.083612040133779306,"#183550"],[0.086956521739130349,"#183551"],[0.090301003344481684,"#183651"],[0.093645484949832727,"#193652"],[0.096989966555184062,"#193652"],[0.10033444816053511,"#193753"],[0.10367892976588644,"#193754"],[0.10702341137123748,"#193854"],[0.11036789297658853,"#1A3855"],[0.11371237458193986,"#1A3955"],[0.11705685618729091,"#1A3956"],[0.12040133779264224,"#1A3956"],[0.12374581939799328,"#1A3A57"],[0.12709030100334462,"#1B3A57"],[0.13043478260869568,"#1B3B58"],[0.13377926421404671,"#1B3B59"],[0.13712374581939804,"#1B3B59"],[0.1404682274247491,"#1C3C5A"],[0.14381270903010043,"#1C3C5A"],[0.14715719063545146,"#1C3D5B"],[0.1505016722408028,"#1C3D5B"],[0.15384615384615385,"#1C3D5C"],[0.15719063545150488,"#1D3E5C"],[0.16053511705685622,"#1D3E5D"],[0.16387959866220728,"#1D3F5D"],[0.16722408026755861,"#1D3F5E"],[0.17056856187290964,"#1D3F5F"],[0.17391304347826098,"#1E405F"],[0.17725752508361203,"#1E4060"],[0.18060200668896306,"#1E4160"],[0.1839464882943144,"#1E4161"],[0.18729096989966545,"#1E4261"],[0.19063545150501679,"#1F4262"],[0.19397993311036782,"#1F4263"],[0.19732441471571915,"#1F4363"],[0.20066889632107021,"#1F4364"],[0.20401337792642155,"#1F4464"],[0.20735785953177258,"#204465"],[0.21070234113712391,"#204465"],[0.21404682274247497,"#204566"],[0.21739130434782603,"#204566"],[0.22073578595317733,"#214667"],[0.22408026755852839,"#214668"],[0.22742474916387972,"#214768"],[0.23076923076923078,"#214769"],[0.23411371237458209,"#214769"],[0.23745819397993315,"#22486A"],[0.2408026755852842,"#22486A"],[0.24414715719063554,"#22496B"],[0.24749163879598657,"#22496C"],[0.2508361204013379,"#224A6C"],[0.25418060200668896,"#234A6D"],[0.2575250836120403,"#234A6D"],[0.26086956521739135,"#234B6E"],[0.26421404682274235,"#234B6E"],[0.26755852842809369,"#244C6F"],[0.27090301003344475,"#244C70"],[0.27424749163879608,"#244C70"],[0.27759197324414714,"#244D71"],[0.28093645484949847,"#244D71"],[0.28428093645484953,"#254E72"],[0.28762541806020053,"#254E72"],[0.29096989966555187,"#254F73"],[0.29431438127090293,"#254F74"],[0.29765886287625426,"#254F74"],[0.30100334448160532,"#265075"],[0.30434782608695665,"#265075"],[0.30769230769230771,"#265176"],[0.31103678929765871,"#265176"],[0.31438127090301005,"#275277"],[0.3177257525083611,"#275278"],[0.32107023411371244,"#275278"],[0.3244147157190635,"#275379"],[0.32775919732441483,"#275379"],[0.33110367892976589,"#28547A"],[0.33444816053511695,"#28547B"],[0.33779264214046828,"#28557B"],[0.34113712374581956,"#28557C"],[0.34448160535117062,"#28567C"],[0.34782608695652167,"#29567D"],[0.35117056856187301,"#29567D"],[0.35451505016722407,"#29577E"],[0.35785953177257512,"#29577F"],[0.36120401337792646,"#2A587F"],[0.36454849498327779,"#2A5880"],[0.3678929765886288,"#2A5980"],[0.37123745819397985,"#2A5981"],[0.37458193979933119,"#2A5982"],[0.37792642140468224,"#2B5A82"],[0.38127090301003358,"#2B5A83"],[0.38461538461538464,"#2B5B83"],[0.38795986622073597,"#2B5B84"],[0.39130434782608697,"#2C5C85"],[0.39464882943143803,"#2C5C85"],[0.39799331103678937,"#2C5D86"],[0.40133779264214042,"#2C5D86"],[0.40468227424749176,"#2C5D87"],[0.40802675585284282,"#2D5E87"],[0.41137123745819415,"#2D5E88"],[0.41471571906354515,"#2D5F89"],[0.41806020066889621,"#2D5F89"],[0.42140468227424754,"#2E608A"],[0.4247491638795986,"#2E608A"],[0.42809364548494994,"#2E618B"],[0.43143812709030099,"#2E618C"],[0.43478260869565233,"#2E618C"],[0.43812709030100333,"#2F628D"],[0.44147157190635439,"#2F628D"],[0.44481605351170572,"#2F638E"],[0.44816053511705678,"#2F638F"],[0.45150501672240811,"#30648F"],[0.45484949832775917,"#306490"],[0.45819397993311051,"#306590"],[0.46153846153846156,"#306591"],[0.46488294314381257,"#306592"],[0.4682274247491639,"#316692"],[0.47157190635451524,"#316693"],[0.47491638795986629,"#316793"],[0.47826086956521735,"#316794"],[0.48160535117056869,"#326895"],[0.48494983277591974,"#326895"],[0.48829431438127074,"#326996"],[0.49163879598662208,"#326996"],[0.49498327759197341,"#326997"],[0.49832775919732447,"#336A98"],[0.50167224080267547,"#336A98"],[0.50501672240802686,"#336B99"],[0.50836120401337792,"#336B99"],[0.5117056856187292,"#346C9A"],[0.51505016722408026,"#346C9B"],[0.51839464882943165,"#346D9B"],[0.52173913043478271,"#346D9C"],[0.52508361204013365,"#346E9D"],[0.52842809364548504,"#356E9D"],[0.5317725752508361,"#356E9E"],[0.53511705685618738,"#356F9E"],[0.53846153846153844,"#356F9F"],[0.54180602006688983,"#3670A0"],[0.54515050167224088,"#3670A0"],[0.54849498327759194,"#3671A1"],[0.55183946488294322,"#3671A1"],[0.55518394648829428,"#3772A2"],[0.55852842809364556,"#3772A3"],[0.56187290969899661,"#3773A3"],[0.565217391304348,"#3773A4"],[0.56856187290969906,"#3773A4"],[0.57190635451505012,"#3874A5"],[0.5752508361204014,"#3874A6"],[0.57859531772575246,"#3875A6"],[0.58193979933110374,"#3875A7"],[0.58528428093645479,"#3976A8"],[0.58862876254180618,"#3976A8"],[0.59197324414715724,"#3977A9"],[0.5953177257525083,"#3977A9"],[0.59866220735785958,"#3978AA"],[0.60200668896321097,"#3A78AB"],[0.60535117056856191,"#3A79AB"],[0.60869565217391297,"#3A79AC"],[0.61204013377926436,"#3A79AC"],[0.61538461538461542,"#3B7AAD"],[0.61872909698996648,"#3B7AAE"],[0.62207357859531776,"#3B7BAE"],[0.62541806020066915,"#3B7BAF"],[0.62876254180602009,"#3C7CB0"],[0.63210702341137115,"#3C7CB0"],[0.63545150501672254,"#3C7DB1"],[0.6387959866220736,"#3C7DB1"],[0.64214046822742465,"#3C7EB2"],[0.64548494983277593,"#3D7EB3"],[0.64882943143812732,"#3D7FB3"],[0.65217391304347827,"#3D7FB4"],[0.65551839464882933,"#3D7FB5"],[0.65886287625418072,"#3E80B5"],[0.66220735785953178,"#3E80B6"],[0.66555183946488283,"#3E81B6"],[0.66889632107023411,"#3E81B7"],[0.6722408026755855,"#3F82B8"],[0.67558528428093656,"#3F82B8"],[0.67892976588628751,"#3F83B9"],[0.6822742474916389,"#3F83BA"],[0.68561872909699018,"#4084BA"],[0.68896321070234101,"#4084BB"],[0.69230769230769229,"#4085BB"],[0.69565217391304368,"#4085BC"],[0.69899665551839474,"#4086BD"],[0.70234113712374568,"#4186BD"],[0.70568561872909707,"#4186BE"],[0.70903010033444835,"#4187BF"],[0.71237458193979974,"#4187BF"],[0.71571906354515047,"#4288C0"],[0.71906354515050186,"#4288C1"],[0.72240802675585258,"#4289C1"],[0.72575250836120386,"#4289C2"],[0.72909698996655525,"#438AC2"],[0.73244147157190653,"#438AC3"],[0.73578595317725792,"#438BC4"],[0.73913043478260865,"#438BC4"],[0.74247491638796004,"#438CC5"],[0.74581939799331076,"#448CC6"],[0.74916387959866204,"#448DC6"],[0.75250836120401343,"#448DC7"],[0.75585284280936471,"#448EC8"],[0.7591973244147161,"#458EC8"],[0.76254180602006683,"#458FC9"],[0.76588628762541822,"#458FC9"],[0.76923076923076894,"#458FCA"],[0.77257525083612033,"#4690CB"],[0.77591973244147161,"#4690CB"],[0.77926421404682289,"#4691CC"],[0.78260869565217428,"#4691CD"],[0.785953177257525,"#4792CD"],[0.78929765886287639,"#4792CE"],[0.79264214046822712,"#4793CF"],[0.79598662207357851,"#4793CF"],[0.79933110367892979,"#4894D0"],[0.80267558528428107,"#4894D0"],[0.80602006688963246,"#4895D1"],[0.80936454849498318,"#4895D2"],[0.81270903010033457,"#4896D2"],[0.81605351170568585,"#4996D3"],[0.81939799331103669,"#4997D4"],[0.82274247491638797,"#4997D4"],[0.82608695652173936,"#4998D5"],[0.82943143812709064,"#4A98D6"],[0.83277591973244136,"#4A99D6"],[0.83612040133779275,"#4A99D7"],[0.83946488294314403,"#4A9AD8"],[0.84280936454849487,"#4B9AD8"],[0.84615384615384615,"#4B9BD9"],[0.84949832775919754,"#4B9BDA"],[0.85284280936454882,"#4B9BDA"],[0.85618729096989954,"#4C9CDB"],[0.85953177257525093,"#4C9CDB"],[0.86287625418060221,"#4C9DDC"],[0.86622073578595304,"#4C9DDD"],[0.86956521739130432,"#4D9EDD"],[0.87290969899665571,"#4D9EDE"],[0.87625418060200699,"#4D9FDF"],[0.87959866220735772,"#4D9FDF"],[0.88294314381270911,"#4DA0E0"],[0.88628762541806039,"#4EA0E1"],[0.88963210702341122,"#4EA1E1"],[0.8929765886287625,"#4EA1E2"],[0.89632107023411389,"#4EA2E3"],[0.89966555183946517,"#4FA2E3"],[0.9030100334448159,"#4FA3E4"],[0.90635451505016729,"#4FA3E5"],[0.90969899665551857,"#4FA4E5"],[0.9130434782608694,"#50A4E6"],[0.91638795986622068,"#50A5E7"],[0.91973244147157207,"#50A5E7"],[0.92307692307692335,"#50A6E8"],[0.92642140468227407,"#51A6E8"],[0.92976588628762546,"#51A7E9"],[0.93311036789297674,"#51A7EA"],[0.93645484949832758,"#51A8EA"],[0.93979933110367886,"#52A8EB"],[0.94314381270903025,"#52A9EC"],[0.94648829431438153,"#52A9EC"],[0.94983277591973225,"#52AAED"],[0.95317725752508364,"#53AAEE"],[0.95652173913043492,"#53ABEE"],[0.95986622073578576,"#53ABEF"],[0.96321070234113704,"#53ACF0"],[0.96655518394648843,"#54ACF0"],[0.96989966555183971,"#54ADF1"],[0.97324414715719043,"#54ADF2"],[0.97658862876254182,"#54AEF2"],[0.9799331103678931,"#55AEF3"],[0.98327759197324394,"#55AFF4"],[0.98662207357859522,"#55AFF4"],[0.98996655518394661,"#55B0F5"],[0.99331103678929789,"#56B0F6"],[0.99665551839464928,"#56B1F6"],[1,"#56B1F7"]],"colorbar":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"thickness":23.039999999999996,"title":"as.numeric(colData(se.subset)$RIN)","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":19.925280199252807},"tickmode":"array","ticktext":["6","7","8"],"tickvals":[0.065967741935483909,0.38747311827956982,0.70897849462365581],"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":17.268576172685766},"ticklen":2,"len":0.5,"yanchor":"top","y":1}},"xaxis":"x","yaxis":"y","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":48.152760481527615,"l":66.749688667496898},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.5043691193667541,0.48612158345155737],"tickmode":"array","ticktext":["-0.50","-0.25","0.00","0.25"],"tickvals":[-0.5,-0.25,0,0.25],"categoryorder":"array","categoryarray":["-0.50","-0.25","0.00","0.25"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":18.596928185969286},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"PC1 (35%)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.53958377357698029,0.59146891804297774],"tickmode":"array","ticktext":["-0.3","0.0","0.3"],"tickvals":[-0.30000000000000004,0,0.30000000000000004],"categoryorder":"array","categoryarray":["-0.3","0.0","0.3"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":18.596928185969286},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"PC2 (14.3%)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"rgba(255,255,255,1)","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":17.268576172685766},"y":0.5,"yanchor":"top","title":{"text":"as.numeric(colData(se.subset)$RIN)<br />Number_Outliers_Binned","font":{"color":"rgba(0,0,0,1)","family":"","size":19.925280199252807}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"14fc2530200c7":{"x":{},"y":{},"colour":{},"label":{},"shape":{},"type":"scatter"}},"cur_data":"14fc2530200c7","visdat":{"14fc2530200c7":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

Top two principal components (PCs) of z-score scaled QC metrics. Shape illustrates number of outliers binned and color shows range of RIN values. 

### QC metrics association to PC components 1: 10


```r
cor1 <- cor(pc$rotation[,1:10], t(qcMetrics.z), method = "pearson")

cor1Res <- apply(cor1, 2, function(x){
  return(names(which.max(abs(x))))
})

uniquePCs <- unique(cor1Res)[order(unique(cor1Res))]
qcpc.list = list()
for(x in 1:length(uniquePCs)){
  pc.num <- uniquePCs[x]
  qcpc.list[[pc.num]] = names(cor1Res)[cor1Res == pc.num]
}

QCPC.table <- plyr::ldply(qcpc.list, cbind)
colnames(QCPC.table) <- c("PC", "Metric")

 QCPC.table %>%
   knitr::kable(
     format = "html", align = rep("c", ncol(QCPC.table) + 1),
     row.names = TRUE
   ) %>%
   kableExtra::kable_styling() %>%
   kableExtra::scroll_box(height = "200px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:200px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> PC </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> Metric </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> 1 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> samtools_reads_properly_paired </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 2 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> rnaseqc_Total.Reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 3 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> rnaseqc_Mapped.Reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 4 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> rnaseqc_Mapped.Unique.Reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 5 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> rnaseqc_High.Quality.Reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 6 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> rnaseqc_rRNA.Reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> rnaseqc_Median.Exon.CV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> STAR_total_reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> STAR_uniquely_mapped </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> STAR_uniquely_mapped_percent </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 11 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> STAR_multimapped_multiple </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> TIN_median </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> rnaseqc_rRNA.Rate </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> rnaseqc_Genes.Detected </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> STAR_avg_mapped_read_length </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> RIN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 17 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> DV200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 18 </td>
   <td style="text-align:center;"> PC3 </td>
   <td style="text-align:center;"> rnaseqc_Unique.Rate.of.Mapped </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 19 </td>
   <td style="text-align:center;"> PC3 </td>
   <td style="text-align:center;"> rnaseqc_Duplicate.Rate.of.Mapped </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20 </td>
   <td style="text-align:center;"> PC3 </td>
   <td style="text-align:center;"> STAR_multimapped_multiple_percent </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 21 </td>
   <td style="text-align:center;"> PC4 </td>
   <td style="text-align:center;"> rnaseqc_Median.of.Avg.Transcript.Coverage </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 22 </td>
   <td style="text-align:center;"> PC4 </td>
   <td style="text-align:center;"> STAR_num_splices </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 23 </td>
   <td style="text-align:center;"> PC5 </td>
   <td style="text-align:center;"> STAR_unmapped_tooshort_percent </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 24 </td>
   <td style="text-align:center;"> PC6 </td>
   <td style="text-align:center;"> rnaseqc_Median.3..bias </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25 </td>
   <td style="text-align:center;"> PC6 </td>
   <td style="text-align:center;"> heterozygosity_mean </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 26 </td>
   <td style="text-align:center;"> PC7 </td>
   <td style="text-align:center;"> rnaseqc_Low.Quality.Reads </td>
  </tr>
</tbody>
</table></div>


Of all of the QC metrics, 60 were most associated, either positively or negatively, to PC1. 


### Heatmap QC outliers

The numeric QC metric values were z-scored across samples, plotted by heatmap. A total of 26 were used. 


```r
Number_Outliers_Binned <- as.factor(qualAnnot)
Number_Outliers_Binned <- factor(Number_Outliers_Binned, 
                                 levels = c("x =< 6",
                                            "6 < x < 20",
                                            "x >= 20"))

annotation_col <- data.frame(Number_Outliers_Binned)

rownames(annotation_col) <- colnames(qcMetrics.z.trim)

annotation_row <- data.frame(Correlated_PC = cor1Res)
rownames(annotation_row) <- rownames(qcMetrics.z.trim)


out <- pheatmap::pheatmap(qcMetrics.z.trim,
  #color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  color = colorRampPalette(c("blue","white","red"), space="rgb")(255),
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  legend = T,
  show_rownames = T,
  show_colnames = T,
  fontsize_row = 7,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "ward.D2"
)
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/heatmap_qc_outiers-1.png" style="display: block; margin: auto;" />

```r
out
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/heatmap_qc_outiers-2.png" style="display: block; margin: auto;" />

 

## Using cutree to assign clusters to samples {.tabset .tabset-fade}
The `cutree` function is used to cut the dendrogram into clusters. Users may define a specific cluster as those with poor quality samples.

### K = 2


```r
##Uses cutree to determine which samples to remove

clusterLabels <- as.data.frame(sort(cutree(out$tree_col, k=2)))
colnames(clusterLabels) <- "Cluster"
clusterLabels %>%
  knitr::kable(
    format = "html", align = rep("c", ncol(QCPC.table) + 1),
    row.names = TRUE
  ) %>%
  kableExtra::kable_styling() %>%
  kableExtra::scroll_box(height = "250px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:250px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> Cluster </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P01_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P02_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P05_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P06_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P09_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P03_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P04_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P07_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P08_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P10_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P11_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P12_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
</tbody>
</table></div>

```r
#  if(output){
#     write.table(clusterLabels, "cutreeClusterLabels_K2.txt", sep = "\t", quote = F)
# }
```



### K = 3


```r
clusterLabels <- as.data.frame(sort(cutree(out$tree_col, k=3)))
colnames(clusterLabels) <- "Cluster"
clusterLabels %>%
  knitr::kable(
    format = "html", align = rep("c", ncol(QCPC.table) + 1),
    row.names = TRUE
  ) %>%
  kableExtra::kable_styling() %>%
  kableExtra::scroll_box(height = "250px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:250px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> Cluster </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P01_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P02_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P05_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P06_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P09_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P03_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P07_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P04_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P08_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P10_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P11_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P12_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
</tbody>
</table></div>

```r
# if(output){
#     write.table(clusterLabels, "cutreeClusterLabels_K3.txt", sep = "\t", quote = F)
# }
```


### K = 4


```r
clusterLabels <- as.data.frame(sort(cutree(out$tree_col, k=4)))
colnames(clusterLabels) <- "Cluster"
clusterLabels %>%
  knitr::kable(
    format = "html", align = rep("c", ncol(QCPC.table) + 1),
    row.names = TRUE
  ) %>%
  kableExtra::kable_styling() %>%
  kableExtra::scroll_box(height = "250px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:250px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> Cluster </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P01_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P02_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P05_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P06_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P09_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P03_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P07_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P04_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P10_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P11_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P12_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P08_S1 </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
</tbody>
</table></div>

```r
# if(output){
#     write.table(clusterLabels, "cutreeClusterLabels_K4.txt", sep = "\t", quote = F)
# }
```


### K = 5


```r
clusterLabels <- as.data.frame(sort(cutree(out$tree_col, k=5)))
colnames(clusterLabels) <- "Cluster"
clusterLabels %>%
  knitr::kable(
    format = "html", align = rep("c", ncol(QCPC.table) + 1),
    row.names = TRUE
  ) %>%
  kableExtra::kable_styling() %>%
  kableExtra::scroll_box(height = "250px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:250px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> Cluster </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P01_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P02_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P05_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P03_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P07_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P04_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P10_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P11_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P12_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P06_S1 </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P09_S1 </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P08_S1 </td>
   <td style="text-align:center;"> 5 </td>
  </tr>
</tbody>
</table></div>

```r
# if(output){
#     write.table(clusterLabels, "cutreeClusterLabels_K5.txt", sep = "\t", quote = F)
# }
```




## QC based sex annotation (based on X and Y-linked genes) {.tabset .tabset-fade}

### Heatmap sex genes 


```r
# -----------------------------
# 0. Build sex-gene matrix
# -----------------------------
control.symbols <- c("XIST", "TXLNGY", "DDX3Y", "KDM5D", "RPS4Y1", "USP9Y", "UTY")

# expression for biopsy subset only
data.dge.sex <- as.matrix(data.cpm)
mode(data.dge.sex) <- "numeric"

gene_symbols <- as.character(mcols(se.subset)$external_gene_name)

sex.df <- data.frame(
  external_gene_name = gene_symbols,
  data.dge.sex,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sex.genes.df <- sex.df %>%
  filter(external_gene_name %in% control.symbols)

sex.expr <- as.matrix(sex.genes.df[, setdiff(colnames(sex.genes.df), "external_gene_name"), drop = FALSE])
mode(sex.expr) <- "numeric"
rownames(sex.expr) <- sex.genes.df$external_gene_name

# collapse duplicate symbols if present
if (any(duplicated(rownames(sex.expr)))) {
  sex.expr <- rowsum(sex.expr, group = rownames(sex.expr), na.rm = TRUE) /
    as.vector(table(rownames(sex.expr)))
}

# keep only genes present
keep_genes <- intersect(control.symbols, rownames(sex.expr))
sex.expr <- sex.expr[keep_genes, , drop = FALSE]

# remove all-NA rows
sex.expr <- sex.expr[rowSums(is.na(sex.expr)) < ncol(sex.expr), , drop = FALSE]

# impute remaining NAs gene-wise by median
if (nrow(sex.expr) > 0) {
  for (i in seq_len(nrow(sex.expr))) {
    na_idx <- is.na(sex.expr[i, ])
    if (any(na_idx)) {
      sex.expr[i, na_idx] <- median(sex.expr[i, ], na.rm = TRUE)
    }
  }
}


sex.expr.z<-t(scale(t(sex.expr),center=T,scale=T))
# BUG FIX from the original script: this used to hardcode
# `c(rep("Male", 6), "Female")`, i.e. it assumed exactly 6 Y-linked genes were
# present, in a fixed order, followed by XIST -- silently mislabeling rows
# whenever a marker gene was absent from the SE (a real risk on any dataset
# other than the original one). Build the row-level marker-type annotation
# from the actual surviving gene rownames instead.
marker_type <- ifelse(rownames(sex.expr.z) %in% "XIST", "X-linked (XIST)", "Y-linked")
annot_row <- data.frame(marker_type = factor(marker_type, levels = c("Y-linked", "X-linked (XIST)")))
rownames(annot_row) <- rownames(sex.expr.z)
sex.expr.z <- t(apply(sex.expr.z, 1, function(x) ifelse(x > 1, 1, ifelse(x < -1, -1, x))))
#sex.genes.trim.nums.z <- t(apply(sex.genes.trim.nums.z, 1, function(x) ifelse(x > 1, 1, ifelse(x < -1, -1, x))))
sex_colors <- list(marker_type = c("Y-linked" = "lightgreen", "X-linked (XIST)" = "purple"))
heatmap.sex <- pheatmap(
  sex.expr.z,
 color = colorRampPalette(c("blue","white","red"), space="rgb")(255),
 annotation_row = annot_row,
 annotation_colors = sex_colors,
 legend = T,
 fontsize = 10,
 show_rownames = T,
 show_colnames = T,
 fontsize_row = 7,
 cluster_rows = TRUE,
 cluster_cols = TRUE,
 clustering_method = "ward.D2"
)
heatmap.sex
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/heatmap_sex_genes-1.png" style="display: block; margin: auto;" />

Heatmap shows sex gene expression for all samples.

### PCA of sex-marker genes 

```r
# -----------------------------
# 1. PCA on samples
# -----------------------------
pca_sex <- prcomp(t(sex.expr), center = TRUE, scale. = TRUE)

y_present <- intersect(c("TXLNGY", "DDX3Y", "KDM5D", "RPS4Y1", "USP9Y", "UTY"), rownames(sex.expr))
x_present <- intersect(c("XIST", "Xist"), rownames(sex.expr))

qc.sex <- data.frame(
  sample = colnames(sex.expr),
  sex_PC1 = pca_sex$x[, 1],
  sex_PC2 = pca_sex$x[, 2],
  XIST = if (length(x_present) > 0) as.numeric(sex.expr[x_present[1], ]) else NA_real_,
  Y_mean = if (length(y_present) > 0) colMeans(sex.expr[y_present, , drop = FALSE], na.rm = TRUE) else NA_real_,
  stringsAsFactors = FALSE
)

# -----------------------------
# 2. Infer sex from PC1 by kmeans
# -----------------------------
set.seed(1)
km <- kmeans(qc.sex$sex_PC1, centers = 2, nstart = 50)

qc.sex$cluster <- factor(km$cluster)

cluster_summary <- qc.sex %>%
  group_by(cluster) %>%
  summarise(
    med_XIST = median(XIST, na.rm = TRUE),
    med_Y = median(Y_mean, na.rm = TRUE),
    .groups = "drop"
  )

female_cluster <- cluster_summary$cluster[which.max(cluster_summary$med_XIST)]
male_cluster   <- cluster_summary$cluster[which.max(cluster_summary$med_Y)]

qc.sex <- qc.sex %>%
  mutate(
    inferred_sex = case_when(
      cluster == female_cluster ~ "Female",
      cluster == male_cluster ~ "Male",
      TRUE ~ "Unknown"
    )
  )

# -----------------------------
# 3. Confidence / outlier metrics
# -----------------------------
sex_medians <- tapply(qc.sex$sex_PC1, qc.sex$inferred_sex, median, na.rm = TRUE)
sex_medians <- sex_medians[names(sex_medians) %in% c("Female", "Male")]

boundary <- mean(sex_medians)
qc.sex$boundary_dist <- abs(qc.sex$sex_PC1 - boundary)

qc.sex$mdist <- NA_real_

for (sx in unique(qc.sex$inferred_sex)) {
  idx <- which(qc.sex$inferred_sex == sx)
  if (length(idx) >= 3) {
    coords <- qc.sex[idx, c("sex_PC1", "sex_PC2"), drop = FALSE]
    center <- colMeans(coords)
    covmat <- cov(coords)

    if (all(is.finite(covmat)) && det(as.matrix(covmat)) != 0) {
      qc.sex$mdist[idx] <- mahalanobis(coords, center = center, cov = covmat)
    }
  }
}

pc1_margin_thresh <- quantile(qc.sex$boundary_dist, 0.10, na.rm = TRUE)
mdist_thresh <- quantile(qc.sex$mdist, 0.95, na.rm = TRUE)

qc.sex <- qc.sex %>%
  mutate(
    flag_sex_near_boundary = boundary_dist <= pc1_margin_thresh,
    flag_sex_within_group_outlier = mdist >= mdist_thresh,
    flag_sex_QC_PC_outlier = flag_sex_near_boundary | flag_sex_within_group_outlier,
    sex_qc_flag = case_when(
      flag_sex_near_boundary & flag_sex_within_group_outlier ~ "Ambiguous + outlier",
      flag_sex_near_boundary ~ "Near boundary",
      flag_sex_within_group_outlier ~ "Within-sex outlier",
      TRUE ~ "Pass"
    )
  )

# optional plot
ggplot(qc.sex, aes(sex_PC1, sex_PC2, color = inferred_sex, shape = sex_qc_flag)) +
  geom_point(size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(sex_qc_flag != "Pass", sample, "")),
    size = 3,
    max.overlaps = 50
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Sex QC by PCA of sex-marker genes",
    subtitle = "Flagged samples are near the sex boundary or outliers within their inferred sex group",
    x = "Sex PC1",
    y = "Sex PC2"
  )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/sex_genes_pca-1.png" style="display: block; margin: auto;" />


## Flagged QC samples {.tabset .tabset-fade}

### Flagged samples table


```r
# Pick the sample universe for this tissue subset
samples <- colnames(se.subset)
## Table 
qc_tbl <- tibble(
  sample_id  = samples,
  patient_id = colData(se.subset)$patient_id %>% as.character()
)

# Flag #1: Somalier cross-patient connections (relatedness >= cutoff)
som_edges <- somalier_pairs %>%
  transmute(a = X.sample_a, b = sample_b, rel = relatedness) %>%
  filter(!is.na(rel), rel >= somalier_cutoff) %>%
  filter(a %in% samples, b %in% samples)

pid <- setNames(qc_tbl$patient_id, qc_tbl$sample_id)
cross_som <- pid[som_edges$a] != pid[som_edges$b]

som_flag_samples <- unique(c(som_edges$a[cross_som], som_edges$b[cross_som]))

qc_tbl <- qc_tbl %>%
  mutate(flag_somalier_cross_patient = sample_id %in% som_flag_samples)
# mismatch
som_flag_pairs <- som_edges[cross_som, ] %>%
  mutate(pid_a = pid[a], pid_b = pid[b])

# Flag #2: HLA genotype cross-patient connections (>= cutoff)
hla_edges <- geno.all.pairs %>%
  transmute(a = sample_a, b = sample_b, rel = percentOFcomparison) %>%
  filter(!is.na(rel), rel >= hla_cutoff) %>%
  filter(a %in% samples, b %in% samples)

cross_hla <- pid[hla_edges$a] != pid[hla_edges$b]
hla_flag_samples <- unique(c(hla_edges$a[cross_hla], hla_edges$b[cross_hla]))

qc_tbl <- qc_tbl %>%
  mutate(flag_hla_cross_patient = sample_id %in% hla_flag_samples)

df_qc <- as.data.frame(qcMetrics.subset)

df_qc$sample_id <- rownames(df_qc)


qc_tbl <- qc_tbl[qc_tbl$sample_id %in% df_qc$sample_id,]
# join QC columns into qc_tbl
qc_tbl <- qc_tbl %>%
  left_join(df_qc, by = "sample_id")


qc_tbl <- qc_tbl %>%
  mutate(
    flag_low_tin_rin = (TIN_median < tin_cutoff & RIN < rin_cutoff),
    flag_3prime_bias = (rnaseqc_Median.3..bias > threeprime_bias_cutoff),
    flag_low_genes_or_high_exon_cv = (rnaseqc_Genes.Detected < gene_detected_cutoff |
                                        rnaseqc_Median.Exon.CV > exon_cv_cutoff),
    flag_high_rrna = (rnaseqc_rRNA.Rate > rrna_rate_cutoff),
    flag_low_dv200 = (DV200 < dv200_cutoff),
    flag_high_hetero = (heterozygosity_mean > hetero_cutoff)
  )

# Flag #4: Expression QC (sample-sample corr outliers + PCA outliers)
# you already computed: colnames(gene.cor.tally)[gene.cor.tally > 10]
expr_corr_outliers <- colnames(gene.cor.tally)[gene.cor.tally > 10]

qc_tbl <- qc_tbl %>%
  mutate(flag_expr_corr_outlier = sample_id %in% expr_corr_outliers)

## Expression PCA outliers
expr_pca_outliers <- colData(se.subset)$sample_id[which(out.color == "yes")]

qc_tbl <- qc_tbl %>%
  mutate(flag_expr_pca_outlier = sample_id %in% expr_pca_outliers)

# Flag #5: QC histogram outliers
highAmountOutlier <- colnames(qc.tally)[which(qc.tally > 6)]
qc_tbl <- qc_tbl %>%
  mutate(flag_many_qc_outliers = sample_id %in% highAmountOutlier)

# Flag #6: QC metrics association to PC1:10
# Get sample PC scores

# convert DFrame safely
#df_qc <- as.data.frame(qcMetrics.subset)

# keep only numeric QC metrics
df_qc_num <- df_qc[, sapply(df_qc, is.numeric)]

# remove columns with zero variance
df_qc_num <- df_qc_num[, apply(df_qc_num, 2, sd, na.rm=TRUE) > 0]

# z-score scale QC metrics (column-wise scaling)
z <- scale(df_qc_num)

# PCA
pc <- prcomp(z, center = TRUE, scale. = FALSE)

#Flag samples ±2 SD
# sample PC scores
pc_scores <- as.data.frame(pc$x[,1:2])
pc_scores$sample_id <- rownames(pc_scores)

# helper
flag_sd <- function(v, k = 2){
  abs(v - mean(v, na.rm = TRUE)) > k * sd(v, na.rm = TRUE)
}

# flags
pc_scores$flag_PC1_2SD <- flag_sd(pc_scores$PC1)
pc_scores$flag_PC2_2SD <- flag_sd(pc_scores$PC2)

# final QC-PCA flag (Flag #7)
pc_scores$flag_QC_PC_outlier <- pc_scores$flag_PC1_2SD | pc_scores$flag_PC2_2SD


qc_tbl$flag_QC_PC_outlier <- pc_scores$flag_QC_PC_outlier[
  match(qc_tbl$sample_id, pc_scores$sample_id)
]

qc_tbl$QC_PC1 <- pc_scores$PC1[match(qc_tbl$sample_id, pc_scores$sample_id)]
qc_tbl$QC_PC2 <- pc_scores$PC2[match(qc_tbl$sample_id, pc_scores$sample_id)]

# Flag #7: Sex QC PCA
# merge sex QC into qc_tbl
qc_tbl <- qc_tbl %>%
  left_join(
    qc.sex %>%
      select(
        sample_id = sample,
        inferred_sex,
        sex_PC1,
        sex_PC2,
        XIST,
        Y_mean,
        boundary_dist,
        mdist,
        sex_qc_flag,
        flag_sex_near_boundary,
        flag_sex_within_group_outlier,
        flag_sex_QC_PC_outlier
      ),
    by = "sample_id"
  )

#### Final: a per-sample tally + per-category counts

flag_cols <- grep("^flag_", names(qc_tbl), value = TRUE)

qc_tbl <- qc_tbl %>%
  mutate(
    n_flags = rowSums(across(all_of(flag_cols), ~ as.integer(replace_na(.x, FALSE))), na.rm = TRUE)
  )

qc_tbl <- qc_tbl %>%
  column_to_rownames("sample_id")

qc_tbl %>%
  knitr::kable(
    format = "html", align = rep("c", ncol(QCPC.table) + 1),
    row.names = TRUE
  ) %>%
  kableExtra::kable_styling() %>%
  kableExtra::scroll_box(height = "250px")
```

<div style="border: 1px solid #ddd; padding: 0px; overflow-y: scroll; height:250px; "><table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;">   </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> patient_id </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_somalier_cross_patient </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_hla_cross_patient </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> samtools_reads_properly_paired </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> TIN_median </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Total.Reads </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Mapped.Reads </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Mapped.Unique.Reads </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Unique.Rate.of.Mapped </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_High.Quality.Reads </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Low.Quality.Reads </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_rRNA.Reads </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_rRNA.Rate </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Median.3..bias </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Genes.Detected </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Median.Exon.CV </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Duplicate.Rate.of.Mapped </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> rnaseqc_Median.of.Avg.Transcript.Coverage </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> STAR_total_reads </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> STAR_uniquely_mapped </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> STAR_uniquely_mapped_percent </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> STAR_multimapped_multiple </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> STAR_multimapped_multiple_percent </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> STAR_unmapped_tooshort_percent </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> STAR_avg_mapped_read_length </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> STAR_num_splices </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> RIN </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> DV200 </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> heterozygosity_mean </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_low_tin_rin </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_3prime_bias </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_low_genes_or_high_exon_cv </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_high_rrna </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_low_dv200 </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_high_hetero </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_expr_corr_outlier </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_expr_pca_outlier </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_many_qc_outliers </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_QC_PC_outlier </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> QC_PC1 </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> QC_PC2 </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> inferred_sex </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> sex_PC1 </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> sex_PC2 </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> XIST </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> Y_mean </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> boundary_dist </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> mdist </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> sex_qc_flag </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_sex_near_boundary </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_sex_within_group_outlier </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> flag_sex_QC_PC_outlier </th>
   <th style="text-align:center;position: sticky; top:0; background-color: #FFFFFF;"> n_flags </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P01_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P01 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 25115640 </td>
   <td style="text-align:center;"> 60.2 </td>
   <td style="text-align:center;"> 30758024 </td>
   <td style="text-align:center;"> 27377041 </td>
   <td style="text-align:center;"> 19481919 </td>
   <td style="text-align:center;"> 0.712 </td>
   <td style="text-align:center;"> 26973631 </td>
   <td style="text-align:center;"> 3784393 </td>
   <td style="text-align:center;"> 575953 </td>
   <td style="text-align:center;"> 0.0187 </td>
   <td style="text-align:center;"> 0.698 </td>
   <td style="text-align:center;"> 11637 </td>
   <td style="text-align:center;"> 0.472 </td>
   <td style="text-align:center;"> 0.289 </td>
   <td style="text-align:center;"> 24.68 </td>
   <td style="text-align:center;"> 30758024 </td>
   <td style="text-align:center;"> 27913156 </td>
   <td style="text-align:center;"> 90.75 </td>
   <td style="text-align:center;"> 544060 </td>
   <td style="text-align:center;"> 1.77 </td>
   <td style="text-align:center;"> 4.90 </td>
   <td style="text-align:center;"> 110.9 </td>
   <td style="text-align:center;"> 4124070 </td>
   <td style="text-align:center;"> 5.8 </td>
   <td style="text-align:center;"> 87.4 </td>
   <td style="text-align:center;"> 0.684 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> -2.3005344 </td>
   <td style="text-align:center;"> 2.2740095 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.492984 </td>
   <td style="text-align:center;"> -0.4517857 </td>
   <td style="text-align:center;"> 14.046702 </td>
   <td style="text-align:center;"> 6.759205 </td>
   <td style="text-align:center;"> 2.488336 </td>
   <td style="text-align:center;"> 3.3437428 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P02_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P02 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 24842441 </td>
   <td style="text-align:center;"> 45.2 </td>
   <td style="text-align:center;"> 27246733 </td>
   <td style="text-align:center;"> 25460582 </td>
   <td style="text-align:center;"> 20482123 </td>
   <td style="text-align:center;"> 0.804 </td>
   <td style="text-align:center;"> 23216815 </td>
   <td style="text-align:center;"> 4029918 </td>
   <td style="text-align:center;"> 270911 </td>
   <td style="text-align:center;"> 0.0099 </td>
   <td style="text-align:center;"> 0.628 </td>
   <td style="text-align:center;"> 12426 </td>
   <td style="text-align:center;"> 0.669 </td>
   <td style="text-align:center;"> 0.144 </td>
   <td style="text-align:center;"> 37.50 </td>
   <td style="text-align:center;"> 27246733 </td>
   <td style="text-align:center;"> 24263011 </td>
   <td style="text-align:center;"> 89.05 </td>
   <td style="text-align:center;"> 552749 </td>
   <td style="text-align:center;"> 2.03 </td>
   <td style="text-align:center;"> 1.54 </td>
   <td style="text-align:center;"> 129.6 </td>
   <td style="text-align:center;"> 1044598 </td>
   <td style="text-align:center;"> 8.0 </td>
   <td style="text-align:center;"> 58.8 </td>
   <td style="text-align:center;"> 0.699 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> -3.3946606 </td>
   <td style="text-align:center;"> -0.1267608 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.578991 </td>
   <td style="text-align:center;"> 0.0544227 </td>
   <td style="text-align:center;"> 5.461023 </td>
   <td style="text-align:center;"> 12.408489 </td>
   <td style="text-align:center;"> 2.583639 </td>
   <td style="text-align:center;"> 1.0779157 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P03_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P03 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 38284044 </td>
   <td style="text-align:center;"> 88.2 </td>
   <td style="text-align:center;"> 40743045 </td>
   <td style="text-align:center;"> 39681127 </td>
   <td style="text-align:center;"> 36248724 </td>
   <td style="text-align:center;"> 0.914 </td>
   <td style="text-align:center;"> 36600299 </td>
   <td style="text-align:center;"> 4142746 </td>
   <td style="text-align:center;"> 507991 </td>
   <td style="text-align:center;"> 0.0125 </td>
   <td style="text-align:center;"> 0.677 </td>
   <td style="text-align:center;"> 9567 </td>
   <td style="text-align:center;"> 1.025 </td>
   <td style="text-align:center;"> 0.233 </td>
   <td style="text-align:center;"> 13.07 </td>
   <td style="text-align:center;"> 40743045 </td>
   <td style="text-align:center;"> 31902029 </td>
   <td style="text-align:center;"> 78.30 </td>
   <td style="text-align:center;"> 702788 </td>
   <td style="text-align:center;"> 1.72 </td>
   <td style="text-align:center;"> 1.58 </td>
   <td style="text-align:center;"> 108.7 </td>
   <td style="text-align:center;"> 4761235 </td>
   <td style="text-align:center;"> 5.9 </td>
   <td style="text-align:center;"> 54.0 </td>
   <td style="text-align:center;"> 0.629 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 1.3937332 </td>
   <td style="text-align:center;"> 2.3224617 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -3.136310 </td>
   <td style="text-align:center;"> 0.4629217 </td>
   <td style="text-align:center;"> 14.288935 </td>
   <td style="text-align:center;"> 5.931384 </td>
   <td style="text-align:center;"> 3.131663 </td>
   <td style="text-align:center;"> 3.5201894 </td>
   <td style="text-align:center;"> Within-sex outlier </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> 6 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P04_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P04 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 38260617 </td>
   <td style="text-align:center;"> 54.7 </td>
   <td style="text-align:center;"> 42511317 </td>
   <td style="text-align:center;"> 41405395 </td>
   <td style="text-align:center;"> 32578082 </td>
   <td style="text-align:center;"> 0.787 </td>
   <td style="text-align:center;"> 40588231 </td>
   <td style="text-align:center;"> 1923086 </td>
   <td style="text-align:center;"> 434245 </td>
   <td style="text-align:center;"> 0.0102 </td>
   <td style="text-align:center;"> 0.606 </td>
   <td style="text-align:center;"> 13822 </td>
   <td style="text-align:center;"> 0.786 </td>
   <td style="text-align:center;"> 0.435 </td>
   <td style="text-align:center;"> 12.76 </td>
   <td style="text-align:center;"> 42511317 </td>
   <td style="text-align:center;"> 32431523 </td>
   <td style="text-align:center;"> 76.29 </td>
   <td style="text-align:center;"> 1236763 </td>
   <td style="text-align:center;"> 2.91 </td>
   <td style="text-align:center;"> 4.09 </td>
   <td style="text-align:center;"> 111.1 </td>
   <td style="text-align:center;"> 4974997 </td>
   <td style="text-align:center;"> 7.5 </td>
   <td style="text-align:center;"> 86.7 </td>
   <td style="text-align:center;"> 0.688 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 0.7502622 </td>
   <td style="text-align:center;"> -0.0361129 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.688555 </td>
   <td style="text-align:center;"> -0.2168243 </td>
   <td style="text-align:center;"> 6.108012 </td>
   <td style="text-align:center;"> 12.633802 </td>
   <td style="text-align:center;"> 2.693203 </td>
   <td style="text-align:center;"> 2.7353319 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P05_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P05 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 21255792 </td>
   <td style="text-align:center;"> 70.7 </td>
   <td style="text-align:center;"> 25166274 </td>
   <td style="text-align:center;"> 22503614 </td>
   <td style="text-align:center;"> 16492012 </td>
   <td style="text-align:center;"> 0.733 </td>
   <td style="text-align:center;"> 20743117 </td>
   <td style="text-align:center;"> 4423157 </td>
   <td style="text-align:center;"> 77193 </td>
   <td style="text-align:center;"> 0.0031 </td>
   <td style="text-align:center;"> 0.591 </td>
   <td style="text-align:center;"> 17702 </td>
   <td style="text-align:center;"> 0.364 </td>
   <td style="text-align:center;"> 0.211 </td>
   <td style="text-align:center;"> 19.71 </td>
   <td style="text-align:center;"> 25166274 </td>
   <td style="text-align:center;"> 22673331 </td>
   <td style="text-align:center;"> 90.09 </td>
   <td style="text-align:center;"> 1027526 </td>
   <td style="text-align:center;"> 4.08 </td>
   <td style="text-align:center;"> 4.24 </td>
   <td style="text-align:center;"> 98.9 </td>
   <td style="text-align:center;"> 2429623 </td>
   <td style="text-align:center;"> 8.6 </td>
   <td style="text-align:center;"> 68.7 </td>
   <td style="text-align:center;"> 0.648 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> -4.4115824 </td>
   <td style="text-align:center;"> -0.7694938 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.360606 </td>
   <td style="text-align:center;"> -0.1245609 </td>
   <td style="text-align:center;"> 13.565740 </td>
   <td style="text-align:center;"> 6.870319 </td>
   <td style="text-align:center;"> 2.355958 </td>
   <td style="text-align:center;"> 0.1620256 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P06_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P06 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 24274042 </td>
   <td style="text-align:center;"> 45.7 </td>
   <td style="text-align:center;"> 30254704 </td>
   <td style="text-align:center;"> 26748838 </td>
   <td style="text-align:center;"> 21228459 </td>
   <td style="text-align:center;"> 0.794 </td>
   <td style="text-align:center;"> 26338478 </td>
   <td style="text-align:center;"> 3916226 </td>
   <td style="text-align:center;"> 172658 </td>
   <td style="text-align:center;"> 0.0057 </td>
   <td style="text-align:center;"> 0.334 </td>
   <td style="text-align:center;"> 17895 </td>
   <td style="text-align:center;"> 1.085 </td>
   <td style="text-align:center;"> 0.335 </td>
   <td style="text-align:center;"> 16.67 </td>
   <td style="text-align:center;"> 30254704 </td>
   <td style="text-align:center;"> 26445092 </td>
   <td style="text-align:center;"> 87.41 </td>
   <td style="text-align:center;"> 336175 </td>
   <td style="text-align:center;"> 1.11 </td>
   <td style="text-align:center;"> 1.01 </td>
   <td style="text-align:center;"> 129.5 </td>
   <td style="text-align:center;"> 3990540 </td>
   <td style="text-align:center;"> 8.7 </td>
   <td style="text-align:center;"> 86.2 </td>
   <td style="text-align:center;"> 0.605 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> -3.0090572 </td>
   <td style="text-align:center;"> -3.1184764 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.469237 </td>
   <td style="text-align:center;"> 0.1860361 </td>
   <td style="text-align:center;"> 6.751394 </td>
   <td style="text-align:center;"> 12.413613 </td>
   <td style="text-align:center;"> 2.473885 </td>
   <td style="text-align:center;"> 1.1222580 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P07_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P07 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 39564156 </td>
   <td style="text-align:center;"> 83.7 </td>
   <td style="text-align:center;"> 48717411 </td>
   <td style="text-align:center;"> 42457663 </td>
   <td style="text-align:center;"> 36422521 </td>
   <td style="text-align:center;"> 0.858 </td>
   <td style="text-align:center;"> 40720807 </td>
   <td style="text-align:center;"> 7996604 </td>
   <td style="text-align:center;"> 510158 </td>
   <td style="text-align:center;"> 0.0105 </td>
   <td style="text-align:center;"> 0.581 </td>
   <td style="text-align:center;"> 9765 </td>
   <td style="text-align:center;"> 0.718 </td>
   <td style="text-align:center;"> 0.435 </td>
   <td style="text-align:center;"> 35.27 </td>
   <td style="text-align:center;"> 48717411 </td>
   <td style="text-align:center;"> 38190326 </td>
   <td style="text-align:center;"> 78.39 </td>
   <td style="text-align:center;"> 1514743 </td>
   <td style="text-align:center;"> 3.11 </td>
   <td style="text-align:center;"> 4.83 </td>
   <td style="text-align:center;"> 101.1 </td>
   <td style="text-align:center;"> 4171636 </td>
   <td style="text-align:center;"> 6.2 </td>
   <td style="text-align:center;"> 53.9 </td>
   <td style="text-align:center;"> 0.720 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2.2543941 </td>
   <td style="text-align:center;"> 3.4499221 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.006522 </td>
   <td style="text-align:center;"> -0.0392765 </td>
   <td style="text-align:center;"> 13.150627 </td>
   <td style="text-align:center;"> 7.313191 </td>
   <td style="text-align:center;"> 2.001875 </td>
   <td style="text-align:center;"> 2.5530350 </td>
   <td style="text-align:center;"> Near boundary </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P08_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P08 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 51470357 </td>
   <td style="text-align:center;"> 56.5 </td>
   <td style="text-align:center;"> 58456397 </td>
   <td style="text-align:center;"> 52136199 </td>
   <td style="text-align:center;"> 41579648 </td>
   <td style="text-align:center;"> 0.798 </td>
   <td style="text-align:center;"> 51025862 </td>
   <td style="text-align:center;"> 7430535 </td>
   <td style="text-align:center;"> 472589 </td>
   <td style="text-align:center;"> 0.0081 </td>
   <td style="text-align:center;"> 0.693 </td>
   <td style="text-align:center;"> 13421 </td>
   <td style="text-align:center;"> 1.281 </td>
   <td style="text-align:center;"> 0.128 </td>
   <td style="text-align:center;"> 11.20 </td>
   <td style="text-align:center;"> 58456397 </td>
   <td style="text-align:center;"> 44569660 </td>
   <td style="text-align:center;"> 76.24 </td>
   <td style="text-align:center;"> 2642975 </td>
   <td style="text-align:center;"> 4.52 </td>
   <td style="text-align:center;"> 1.16 </td>
   <td style="text-align:center;"> 147.3 </td>
   <td style="text-align:center;"> 3823436 </td>
   <td style="text-align:center;"> 6.1 </td>
   <td style="text-align:center;"> 80.6 </td>
   <td style="text-align:center;"> 0.635 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 4.5940822 </td>
   <td style="text-align:center;"> -1.6419847 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.408857 </td>
   <td style="text-align:center;"> -0.1616132 </td>
   <td style="text-align:center;"> 7.740877 </td>
   <td style="text-align:center;"> 12.478722 </td>
   <td style="text-align:center;"> 2.413505 </td>
   <td style="text-align:center;"> 3.0142772 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P09_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P09 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 21033829 </td>
   <td style="text-align:center;"> 80.7 </td>
   <td style="text-align:center;"> 24005634 </td>
   <td style="text-align:center;"> 21996753 </td>
   <td style="text-align:center;"> 19190120 </td>
   <td style="text-align:center;"> 0.872 </td>
   <td style="text-align:center;"> 19746009 </td>
   <td style="text-align:center;"> 4259625 </td>
   <td style="text-align:center;"> 450325 </td>
   <td style="text-align:center;"> 0.0188 </td>
   <td style="text-align:center;"> 0.347 </td>
   <td style="text-align:center;"> 11843 </td>
   <td style="text-align:center;"> 0.583 </td>
   <td style="text-align:center;"> 0.381 </td>
   <td style="text-align:center;"> 22.27 </td>
   <td style="text-align:center;"> 24005634 </td>
   <td style="text-align:center;"> 18527689 </td>
   <td style="text-align:center;"> 77.18 </td>
   <td style="text-align:center;"> 598281 </td>
   <td style="text-align:center;"> 2.49 </td>
   <td style="text-align:center;"> 1.15 </td>
   <td style="text-align:center;"> 143.9 </td>
   <td style="text-align:center;"> 2903300 </td>
   <td style="text-align:center;"> 8.9 </td>
   <td style="text-align:center;"> 79.1 </td>
   <td style="text-align:center;"> 0.659 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> -2.9910640 </td>
   <td style="text-align:center;"> 1.0859968 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.451598 </td>
   <td style="text-align:center;"> -0.0653546 </td>
   <td style="text-align:center;"> 13.335281 </td>
   <td style="text-align:center;"> 6.720068 </td>
   <td style="text-align:center;"> 2.446951 </td>
   <td style="text-align:center;"> 0.0295454 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P10_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P10 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 47049104 </td>
   <td style="text-align:center;"> 53.9 </td>
   <td style="text-align:center;"> 50528908 </td>
   <td style="text-align:center;"> 49018897 </td>
   <td style="text-align:center;"> 42761799 </td>
   <td style="text-align:center;"> 0.872 </td>
   <td style="text-align:center;"> 44375315 </td>
   <td style="text-align:center;"> 6153593 </td>
   <td style="text-align:center;"> 553581 </td>
   <td style="text-align:center;"> 0.0110 </td>
   <td style="text-align:center;"> 0.400 </td>
   <td style="text-align:center;"> 14762 </td>
   <td style="text-align:center;"> 1.148 </td>
   <td style="text-align:center;"> 0.380 </td>
   <td style="text-align:center;"> 20.04 </td>
   <td style="text-align:center;"> 50528908 </td>
   <td style="text-align:center;"> 41754223 </td>
   <td style="text-align:center;"> 82.63 </td>
   <td style="text-align:center;"> 602222 </td>
   <td style="text-align:center;"> 1.19 </td>
   <td style="text-align:center;"> 4.66 </td>
   <td style="text-align:center;"> 146.6 </td>
   <td style="text-align:center;"> 2978618 </td>
   <td style="text-align:center;"> 7.5 </td>
   <td style="text-align:center;"> 89.4 </td>
   <td style="text-align:center;"> 0.624 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2.3006723 </td>
   <td style="text-align:center;"> -1.4848140 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.338610 </td>
   <td style="text-align:center;"> 0.2287625 </td>
   <td style="text-align:center;"> 7.160645 </td>
   <td style="text-align:center;"> 12.296752 </td>
   <td style="text-align:center;"> 2.343257 </td>
   <td style="text-align:center;"> 1.8520239 </td>
   <td style="text-align:center;"> Near boundary </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> 5 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P11_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P11 </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 49263052 </td>
   <td style="text-align:center;"> 51.4 </td>
   <td style="text-align:center;"> 57918654 </td>
   <td style="text-align:center;"> 53078089 </td>
   <td style="text-align:center;"> 44517928 </td>
   <td style="text-align:center;"> 0.839 </td>
   <td style="text-align:center;"> 55619440 </td>
   <td style="text-align:center;"> 2299214 </td>
   <td style="text-align:center;"> 406922 </td>
   <td style="text-align:center;"> 0.0070 </td>
   <td style="text-align:center;"> 0.606 </td>
   <td style="text-align:center;"> 10693 </td>
   <td style="text-align:center;"> 0.382 </td>
   <td style="text-align:center;"> 0.286 </td>
   <td style="text-align:center;"> 24.75 </td>
   <td style="text-align:center;"> 57918654 </td>
   <td style="text-align:center;"> 45400243 </td>
   <td style="text-align:center;"> 78.39 </td>
   <td style="text-align:center;"> 900353 </td>
   <td style="text-align:center;"> 1.55 </td>
   <td style="text-align:center;"> 2.78 </td>
   <td style="text-align:center;"> 133.4 </td>
   <td style="text-align:center;"> 2232210 </td>
   <td style="text-align:center;"> 8.1 </td>
   <td style="text-align:center;"> 70.7 </td>
   <td style="text-align:center;"> 0.646 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2.7027916 </td>
   <td style="text-align:center;"> -1.0650787 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.491446 </td>
   <td style="text-align:center;"> 0.1379603 </td>
   <td style="text-align:center;"> 13.176364 </td>
   <td style="text-align:center;"> 6.656604 </td>
   <td style="text-align:center;"> 2.486798 </td>
   <td style="text-align:center;"> 0.3914619 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Data2Cohort_P12_S1 </td>
   <td style="text-align:center;"> Data2Cohort_P12 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 43641521 </td>
   <td style="text-align:center;"> 66.9 </td>
   <td style="text-align:center;"> 52745388 </td>
   <td style="text-align:center;"> 46600063 </td>
   <td style="text-align:center;"> 37625175 </td>
   <td style="text-align:center;"> 0.807 </td>
   <td style="text-align:center;"> 49028714 </td>
   <td style="text-align:center;"> 3716674 </td>
   <td style="text-align:center;"> 331312 </td>
   <td style="text-align:center;"> 0.0063 </td>
   <td style="text-align:center;"> 0.239 </td>
   <td style="text-align:center;"> 12693 </td>
   <td style="text-align:center;"> 1.186 </td>
   <td style="text-align:center;"> 0.275 </td>
   <td style="text-align:center;"> 27.97 </td>
   <td style="text-align:center;"> 52745388 </td>
   <td style="text-align:center;"> 42709550 </td>
   <td style="text-align:center;"> 80.97 </td>
   <td style="text-align:center;"> 1205743 </td>
   <td style="text-align:center;"> 2.29 </td>
   <td style="text-align:center;"> 1.20 </td>
   <td style="text-align:center;"> 112.2 </td>
   <td style="text-align:center;"> 3780049 </td>
   <td style="text-align:center;"> 7.3 </td>
   <td style="text-align:center;"> 78.5 </td>
   <td style="text-align:center;"> 0.768 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2.1109631 </td>
   <td style="text-align:center;"> -0.8896688 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.455217 </td>
   <td style="text-align:center;"> -0.0106879 </td>
   <td style="text-align:center;"> 6.355858 </td>
   <td style="text-align:center;"> 12.366763 </td>
   <td style="text-align:center;"> 2.459864 </td>
   <td style="text-align:center;"> 0.1981933 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
</tbody>
</table></div>

```r
qc_tbl %>%
  tibble::rownames_to_column("sample_id") %>%
  readr::write_tsv(qcFile)
```

### Heatmap of all QC flags (tissue subset)

```r
# Build the flag matrix
flag_cols <- grep("^flag_", names(qc_tbl), value = TRUE)

flag_matrix <- qc_tbl %>%
  select(all_of(flag_cols)) %>%
  mutate(across(everything(), ~ as.integer(replace_na(.x, FALSE)))) %>%
  as.matrix()

# transpose: flags = rows, samples = columns
flag_matrix <- t(flag_matrix)

# assign sample IDs as column names
colnames(flag_matrix) <- rownames(qc_tbl)

# clean flag labels
rownames(flag_matrix) <- gsub("^flag_", "", rownames(flag_matrix))
rownames(flag_matrix) <- gsub("_", " ", rownames(flag_matrix))

# Order samples by QC severity
sample_order <- qc_tbl %>%
  arrange(desc(n_flags)) %>%
  rownames()

flag_matrix <- flag_matrix[, sample_order, drop = FALSE]

# top annotation: continuous covariates
ha_col <- HeatmapAnnotation(
  n_flags = qc_tbl[sample_order, "n_flags"],
  RIN     = qc_tbl[sample_order, "RIN"],
  TIN     = qc_tbl[sample_order, "TIN_median"],
  DV200   = qc_tbl[sample_order, "DV200"],
  which = "column"
)

# top barplot annotation
ha_bar <- HeatmapAnnotation(
  Flags = anno_barplot(
    qc_tbl[sample_order, "n_flags"],
    gp = gpar(fill = "black"),
    border = FALSE,
    height = unit(2.2, "cm")
  ),
  which = "column"
)

# Group QC flags
qc_groups <- c(
  "somalier cross patient"      = "Identity",
  "hla cross patient"           = "Identity",
  "sex near boundary"           = "Identity",
  "sex within group outlier"    = "Identity",
  "sex QC PC outlier"           = "Identity",
  "high hetero"                 = "Identity",

  "low tin rin"                 = "RNA integrity",
  "3prime bias"                 = "RNA integrity",
  "high rrna"                   = "RNA integrity",
  "low dv200"                   = "RNA integrity",

  "low genes or high exon cv"   = "Library quality",

  "expr corr outlier"           = "Expression outlier",
  "expr pca outlier"            = "Expression outlier",
  "QC PC outlier"               = "Expression outlier",

  "many qc outliers"            = "Global failure"
)

row_split <- qc_groups[rownames(flag_matrix)]

ha_top <- c(ha_col, ha_bar)

flag_col_fun <- colorRamp2(c(0, 1), c("white", "firebrick"))

ht <- Heatmap(
  flag_matrix,
  name = "Flag",
  col = flag_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  top_annotation = ha_top,
  row_split = row_split,
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "RNA-seq QC Flag Matrix",
  rect_gp = gpar(col = "grey85")
)

draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/data2_style/output/example_QCReport_data2style_files/figure-html/flagged_all_samples-1.png" style="display: block; margin: auto;" />

### Heatmap of a priority sample subset (optional)

This section is entirely **optional** and only runs if you supply
`params$priority_list_file` — a spreadsheet (.xlsx) with a priority column
(`params$priority_column`, default `"Priority"`) and one or more ID columns
(`params$priority_id_cols`) that may match `rnaAnnot$sample_id`. It highlights the
same QC flag matrix as above, restricted to samples flagged
`params$priority_value` (default `"High"`) in that sheet — useful if your cohort has
a subset of samples (e.g. a specific site, or samples earmarked for a follow-on
analysis) you want to QC-review first. If you don't have such a list, skip this
section (leave `priority_list_file = NULL`, the default).


```r
priorityIds <- readxl::read_xlsx(priority_list_file)

priorityHigh <- priorityIds[!is.na(priorityIds[[priority_column]]) &
                               priorityIds[[priority_column]] == priority_value, ]
target_ids <- rnaAnnot$sample_id

# Try each configured ID column in turn until one matches a known sample_id
match_priority_id <- function(row) {
  for (col in priority_id_cols) {
    if (col %in% names(row) && !is.na(row[[col]]) && row[[col]] %in% target_ids) {
      return(row[[col]])
    }
  }
  NA_character_
}
priorityHigh$Matched_ID <- apply(priorityHigh, 1, match_priority_id)

highPrio <- qc_tbl[rownames(qc_tbl) %in% priorityHigh$Matched_ID, , drop = FALSE]

flag_cols <- grep("^flag_", names(highPrio), value = TRUE)

flag_matrix <- highPrio %>%
  select(all_of(flag_cols)) %>%
  mutate(across(everything(), ~ as.integer(replace_na(.x, FALSE)))) %>%
  as.matrix()

flag_matrix <- t(flag_matrix)
colnames(flag_matrix) <- rownames(highPrio)

rownames(flag_matrix) <- gsub("^flag_", "", rownames(flag_matrix))
rownames(flag_matrix) <- gsub("_", " ", rownames(flag_matrix))

sample_order <- highPrio %>%
  arrange(desc(n_flags)) %>%
  rownames()

flag_matrix <- flag_matrix[, sample_order, drop = FALSE]

ha_col <- HeatmapAnnotation(
  n_flags = highPrio[sample_order, "n_flags"],
  RIN     = highPrio[sample_order, "RIN"],
  TIN     = highPrio[sample_order, "TIN_median"],
  DV200   = highPrio[sample_order, "DV200"],
  which = "column"
)

ha_bar <- HeatmapAnnotation(
  Flags = anno_barplot(
    highPrio[sample_order, "n_flags"],
    gp = gpar(fill = "black"),
    border = FALSE,
    height = unit(2.2, "cm")
  ),
  which = "column"
)

row_split <- qc_groups[rownames(flag_matrix)]
ha_top <- c(ha_col, ha_bar)

flag_col_fun <- colorRamp2(c(0, 1), c("white", "firebrick"))

ht_pri <- Heatmap(
  flag_matrix,
  name = "Flag",
  col = flag_col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 11, fontface = "bold"),
  top_annotation = ha_top,
  row_split = row_split,
  row_title_gp = gpar(fontsize = 11, fontface = "bold"),
  column_title = "RNA-seq QC Flag Matrix (priority subset)",
  rect_gp = gpar(col = "grey85")
)

draw(ht_pri, heatmap_legend_side = "right", annotation_legend_side = "right")
```



## Save filtered SE with appended RNA annotation
The full SE object with appended RNA annotation will be saved to: annotated_Gene_Expression.rds.

```r
# NOTE: the original script had this saveRDS() call commented out, so `output_se`
# was documented but silently never written. Restored to actually respect the param
# (matching the behavior of the per-tissue save below), guarded so it's a no-op
# unless you explicitly set params$output_se.
if (!is.null(output_se)) {
  saveRDS(se, output_se)
}
```

## Save filtered SE per tissue type
The filtered per tissue type SE object with appended RNA annotation will be saved to the file: all_samples_Gene_Expression.rds.

```r
if (!is.null(output_se_per)){
  saveRDS(se.subset, output_se_per)
}
```

## Session Information

```r
rm(list = ls())
gc()
```

```
##            used  (Mb) gc trigger (Mb) limit (Mb) max used (Mb)
## Ncells  8390710 448.2   16140314  862         NA 16140314  862
## Vcells 14928778 113.9   26728521  204      16384 26728521  204
```

```r
sessionInfo()
```

```
## R version 4.3.0 (2023-04-21)
## Platform: aarch64-apple-darwin20 (64-bit)
## Running under: macOS 15.7.7
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRblas.0.dylib 
## LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## time zone: America/New_York
## tzcode source: internal
## 
## attached base packages:
## [1] grid      stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] htmlwidgets_1.6.4           visNetwork_2.1.4           
##  [3] kableExtra_1.4.0            igraph_2.1.4               
##  [5] lubridate_1.9.2             forcats_1.0.0              
##  [7] stringr_1.5.0               dplyr_1.2.1                
##  [9] purrr_1.0.1                 readr_2.1.4                
## [11] tidyr_1.3.2                 tibble_3.2.1               
## [13] tidyverse_2.0.0             plotly_4.12.0              
## [15] ggrepel_0.9.8               ComplexHeatmap_2.18.0      
## [17] ggpubr_0.6.1                readxl_1.4.2               
## [19] pheatmap_1.0.13             ggplot2_4.0.3              
## [21] circlize_0.4.18             RColorBrewer_1.1-3         
## [23] edgeR_4.0.16                limma_3.58.1               
## [25] rmarkdown_2.22              SummarizedExperiment_1.32.0
## [27] Biobase_2.62.0              GenomicRanges_1.54.1       
## [29] GenomeInfoDb_1.38.8         IRanges_2.36.0             
## [31] S4Vectors_0.40.2            BiocGenerics_0.48.1        
## [33] MatrixGenerics_1.14.0       matrixStats_1.5.0          
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-9            rlang_1.3.0             magrittr_2.0.3         
##  [4] clue_0.3-68             GetoptLong_1.1.1        compiler_4.3.0         
##  [7] mgcv_1.8-42             systemfonts_1.3.2       png_0.1-8              
## [10] vctrs_0.7.3             pkgconfig_2.0.3         shape_1.4.6.1          
## [13] crayon_1.5.2            fastmap_1.1.1           backports_1.4.1        
## [16] XVector_0.42.0          labeling_0.4.2          utf8_1.2.3             
## [19] tzdb_0.4.0              bit_4.0.5               xfun_0.39              
## [22] zlibbioc_1.48.2         cachem_1.0.8            jsonlite_1.8.5         
## [25] highr_0.10              DelayedArray_0.28.0     broom_1.0.4            
## [28] parallel_4.3.0          cluster_2.1.4           R6_2.5.1               
## [31] bslib_0.4.2             stringi_1.7.12          car_3.1-3              
## [34] jquerylib_0.1.4         cellranger_1.1.0        Rcpp_1.1.2             
## [37] iterators_1.0.14        knitr_1.43              splines_4.3.0          
## [40] timechange_0.2.0        Matrix_1.5-4            tidyselect_1.2.1       
## [43] rstudioapi_0.14         abind_1.4-8             yaml_2.3.7             
## [46] doParallel_1.0.17       codetools_0.2-19        plyr_1.8.9             
## [49] lattice_0.21-8          withr_3.0.3             S7_0.2.2               
## [52] evaluate_0.21           xml2_1.3.4              pillar_1.9.0           
## [55] carData_3.0-6           foreach_1.5.2           generics_0.1.3         
## [58] vroom_1.6.3             RCurl_1.98-1.19         hms_1.1.3              
## [61] scales_1.4.0            glue_1.6.2              lazyeval_0.2.3         
## [64] tools_4.3.0             data.table_1.14.8       locfit_1.5-9.12        
## [67] ggsignif_0.6.4          crosstalk_1.2.2         colorspace_2.1-0       
## [70] nlme_3.1-162            GenomeInfoDbData_1.2.11 Formula_1.2-5          
## [73] cli_3.6.6               textshaping_0.3.6       fansi_1.0.4            
## [76] viridisLite_0.4.2       S4Arrays_1.2.1          svglite_2.2.1          
## [79] gtable_0.3.6            rstatix_0.7.2           sass_0.4.6             
## [82] digest_0.6.31           SparseArray_1.2.4       rjson_0.2.23           
## [85] farver_2.1.1            htmltools_0.5.9         lifecycle_1.0.5        
## [88] httr_1.4.6              GlobalOptions_0.1.4     statmod_1.5.0          
## [91] bit64_4.0.5
```

