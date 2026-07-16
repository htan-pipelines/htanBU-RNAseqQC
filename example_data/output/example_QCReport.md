---
title: "htanBU-RNAseqQC: Bulk RNA-seq QC Report"
date: "2026-07-15"
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
  # Only `sample_id` is strictly required to match colnames(se); the rest are used for
  # grouping/plotting and can be left as-is (or NA) if you don't have that field.
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
  tissue: NULL                     # value of `tissue_type` to subset on for section 2 (NULL = use all samples)
  site: NULL                       # value of `collection_site` to subset on for section 2 (NULL = don't filter by site)

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
  output_se_per: NULL      # path to save the tissue/site-subset SE (skipped if NULL)
  qcFile: "qc_flag_summary.tsv"   # path to write the per-sample QC flag table (TSV)
  showSession: TRUE
---

# Synthetic Example QC

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
(rendered to a standalone `somalier.html`, not embedded in this report).

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
```

```
## Warning: package 'edgeR' was built under R version 4.3.2
```

```
## Warning: package 'limma' was built under R version 4.3.1
```

```
## Warning: package 'pheatmap' was built under R version 4.3.3
```

```
## Warning: package 'ggpubr' was built under R version 4.3.3
```

```
## Warning: package 'ComplexHeatmap' was built under R version 4.3.1
```

```
## Warning: package 'igraph' was built under R version 4.3.3
```

```
## Warning: package 'kableExtra' was built under R version 4.3.1
```

```r
# ---- Read in params ----
se <- params$se
rnaAnnot <- params$rnaAnnot
somalier_pairs <- params$somalier_pairs
genotypes <- params$genotypes
column_map <- params$column_map
tissueType <- params$tissueType
tissue <- params$tissue
site <- params$site
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

# Re-order rnaAnnot to match the sample order in the SummarizedExperiment
rnaAnnot <- rnaAnnot[match(colnames(se), rnaAnnot$sample_id), ]
rownames(rnaAnnot) <- gsub(" ", "", rnaAnnot$sample_id)
```


## General Statistics Plots {.tabset .tabset-fade}

### Tissue type


```r
colors <- mk_disc_cols(rnaAnnot$tissue_type, "Set1")

ggplot(rnaAnnot, aes(tissue_type, fill = tissue_type)) +
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
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/tissue_type-1.png" style="display: block; margin: auto;" />

### Tissue type per participant


```r
# tissue per patient ID 
summary.samples <-  rnaAnnot %>% select(sample_id,patient_id,tissue_type, collection_site, cohort) %>%
  dplyr::group_by(patient_id, tissue_type, collection_site, cohort) %>%
  dplyr::count(patient_id)
#summary.samples

ggplot(summary.samples, aes(x=patient_id, y=n))+
  geom_bar(aes(fill=tissue_type),stat="identity",position="stack")+
  labs(y="Total tissue count", x = "Participant")+
  #scale_y_continuous(breaks = round(seq(min(summary.samples$n),max(summary.samples$n), by = 1),1))+
  scale_y_continuous(breaks = round(seq(min(summary.samples$n),15, by = 1),1))+
  scale_fill_manual(values = colors) +
  #theme(axis.text.x = element_text(angle=-90, hjust=.1))
  guides(fill = guide_legend(title = "Tissue type"))+
theme_bw() %+replace%
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 19, family = "serif"),
          axis.ticks = element_blank(),
          #axis.title.y = element_blank(),
          #axis.title.x = element_blank(),
          #axis.text.x = element_blank(),
          axis.text.x = element_text(angle=-90, hjust=.1),
          axis.text.y = element_text(color = "black"),
          legend.position = "top",
          legend.text = element_text(size = 20),
          legend.key.size = unit(1, "char")
    )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/tissue_per_participant-1.png" style="display: block; margin: auto;" />

### Tissue typer per site


```r
colors <- mk_disc_cols(rnaAnnot$collection_site, "Set1")
ggplot(summary.samples, aes(fill=collection_site, y=n, x=tissue_type))+
  geom_bar(stat="identity",position="stack")+
  labs(y="Total tissue count", x = NULL)+
  #scale_y_continuous(breaks = round(seq(min(summary.samples$n),max(summary.samples$n), by = 1),1))+
  #scale_y_continuous(breaks = round(seq(min(summary.samples$n),15, by = 1),1))+
  scale_fill_manual(values = colors) +
  #theme(axis.text.x = element_text(angle=-90, hjust=.1))
  guides(fill = guide_legend(title = "Collection site"))+
theme_bw() %+replace%
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          text = element_text(size = 19, family = "serif"),
          axis.ticks = element_blank(),
          #axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          #axis.text.x = element_blank(),
          #axis.text.x = element_text(angle=-90, hjust=.1),
          axis.text.y = element_text(color = "black"),
          legend.position = "top",
          legend.text = element_text(size = 20),
          legend.key.size = unit(1, "char")
    )
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/tissue_per_site-1.png" style="display: block; margin: auto;" />

### Tissue type per cohort


```r
colors <- mk_disc_cols(rnaAnnot$cohort, "Set1")
ggplot(summary.samples, aes(fill=cohort, y=n, x=tissue_type))+
geom_bar(stat="identity",position="stack")+
labs(y="Total tissue count", x = NULL)+
scale_fill_manual(values = colors) +
guides(fill = guide_legend(title = "Cohort"))+
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
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/tissue_per_cohort-1.png" style="display: block; margin: auto;" />


```r
# All cutoffs come from params$qc_cutoffs (see YAML header) so they can be
# overridden per cohort/tissue without editing this file.
# NOTE: this chunk used to live just above "QC Metric Overview" (after tissue/site
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

# mk_disc_cols() is defined once in the setup chunk and reused here
colorList <- list(
  sample_type = mk_disc_cols(rnaAnnot$tissue_type, "Set1"),
  site        = mk_disc_cols(rnaAnnot$collection_site, "Set1"),
  cohort      = mk_disc_cols(rnaAnnot$cohort, "Set1")
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


ha = HeatmapAnnotation(
  df = data.frame(participant = rnaAnnot$patient_id, sample_type= rnaAnnot$tissue_type,
                  site = rnaAnnot$collection_site, cohort= rnaAnnot$cohort,
                  star_uniquely_mapped_per = qc.metrics$STAR_uniquely_mapped_percent,
                  rnaseqc_median_transcript_coverage = qc.metrics$rnaseqc_Median.of.Avg.Transcript.Coverage,
                  rnaseqc_rRNA_rate = qc.metrics$rnaseqc_rRNA.Rate,
                  rnaseqc_high_quality_reads = qc.metrics$rnaseqc_High.Quality.Reads,
                  samtools_properly_paired = qc.metrics$samtools_reads_properly_paired,
                  ranseqc_low_quality_reads = qc.metrics$rnaseqc_Low.Quality.Reads
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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/somalier_semi_supervised_heatmap-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/somalier_unsupervised_heatmap-1.png" style="display: block; margin: auto;" />

```r
# dev.off()
```

### Somalier relatedness network graph

The interactive sample-relatedness network (colored by patient ID and by TIN) is
generated as a **standalone HTML file** rather than embedded in this report — run
`generate_somalier_network.R` (in this repo) against your `somalier.pairs.tsv` to
produce `somalier.html`. It uses the same `somalier_pairs` relatedness cutoff
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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/HLA_genotype_check_alt-1.png" style="display: block; margin: auto;" />

```r
legend("topright",
       legend = patient_ids_uni,
       fill = colors,
       title = "Patient",
       cex = 0.8)
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/HLA_genotype_check_alt-2.png" style="display: block; margin: auto;" />


# QC: (``SiteA Lung``) Analysis
## QC Metric Overview {.tabset .tabset-fade}


```r
se.subset <- se[,which(colData(se)$tissue_type == tissue)]
se.subset <- se.subset[,which(colData(se.subset)$collection_site == site)]
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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/TIN_vs_RIN-1.png" style="display: block; margin: auto;" />
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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/TIN_vs_Median_3_bias-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/Genes_detected_vs_Median_Exon_CV-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/RIN_vs_rRNA_rate-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/TIN_vs_DV200-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/RIN_vs_DV200-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/Somalier_heterozygosity_vs_DV200-1.png" style="display: block; margin: auto;" />




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


mat_col<-data.frame(sample=as.factor(colData(se.subset)$tissue_type))
rownames(mat_col) <- colnames(se.subset)


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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/qc_metrics_heatmap-1.png" style="display: block; margin: auto;" />


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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/Sample_sample_correlation_heatmap_RIN-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/Sample_sample_correlation_heatmap_TIN-1.png" style="display: block; margin: auto;" />

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
#ggplot(as.data.frame(pca$rotation[,c(1,2)]), aes(PC1,PC2, color=colData(se)$TIN_median, shape=colData(se)$tissue_type)) +
sample_id <- rownames(pca$rotation)
plot2 <- ggplot(as.data.frame(pca$rotation[,c(1,2)]), 
                aes(PC1,PC2,color=colData(se.subset)$TIN_median,label = sample_id)) +
  geom_point(size=3,aes(shape=colData(se.subset)$collection_site)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
#plot2
ggplotly(plot2)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-7467675db6975ab6d056" style="width:1056px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-7467675db6975ab6d056">{"x":{"data":[{"x":[-0.41111781717526369,0.3406716837676157,-0.22213877826446352,0.24694489205281545,-0.39487515083993557,0.32717429264698528,-0.23990188662012524,0.25188371136087773,-0.27756807602035954,0.37892712909185017],"y":[0.40943733013196965,0.11972671630713563,0.20647913029965592,0.19226077096856722,-0.28166110669463779,-0.049867851015437524,-0.61528392684004041,-0.42020469336764826,0.26981755589983103,0.16929607431060212],"text":["PC1: -0.4111178<br />PC2:  0.40943733<br />colData(se.subset)$TIN_median: 63.2<br />sample_id: SynCohort_P01_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.3406717<br />PC2:  0.11972672<br />colData(se.subset)$TIN_median: 83.5<br />sample_id: SynCohort_P02_S1<br />colData(se.subset)$collection_site: SiteA","PC1: -0.2221388<br />PC2:  0.20647913<br />colData(se.subset)$TIN_median: 78.4<br />sample_id: SynCohort_P03_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.2469449<br />PC2:  0.19226077<br />colData(se.subset)$TIN_median: 82.2<br />sample_id: SynCohort_P04_S1<br />colData(se.subset)$collection_site: SiteA","PC1: -0.3948752<br />PC2: -0.28166111<br />colData(se.subset)$TIN_median: 60.2<br />sample_id: SynCohort_P05_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.3271743<br />PC2: -0.04986785<br />colData(se.subset)$TIN_median: 75.6<br />sample_id: SynCohort_P06_S1<br />colData(se.subset)$collection_site: SiteA","PC1: -0.2399019<br />PC2: -0.61528393<br />colData(se.subset)$TIN_median: 40.2<br />sample_id: SynCohort_P07_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.2518837<br />PC2: -0.42020469<br />colData(se.subset)$TIN_median: 83.3<br />sample_id: SynCohort_P08_S1<br />colData(se.subset)$collection_site: SiteA","PC1: -0.2775681<br />PC2:  0.26981756<br />colData(se.subset)$TIN_median: 40.4<br />sample_id: SynCohort_P09_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.3789271<br />PC2:  0.16929607<br />colData(se.subset)$TIN_median: 50.8<br />sample_id: SynCohort_P10_S1<br />colData(se.subset)$collection_site: SiteA"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(53,110,158,1)","rgba(86,177,247,1)","rgba(77,160,224,1)","rgba(84,173,241,1)","rgba(48,101,145,1)","rgba(73,150,211,1)","rgba(19,43,67,1)","rgba(86,176,246,1)","rgba(19,44,68,1)","rgba(34,73,107,1)"],"opacity":1,"size":11.338582677165356,"symbol":"circle","line":{"width":1.8897637795275593,"color":["rgba(53,110,158,1)","rgba(86,177,247,1)","rgba(77,160,224,1)","rgba(84,173,241,1)","rgba(48,101,145,1)","rgba(73,150,211,1)","rgba(19,43,67,1)","rgba(86,176,246,1)","rgba(19,44,68,1)","rgba(34,73,107,1)"]}},"hoveron":"points","name":"SiteA","legendgroup":"SiteA","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.40000000000000002],"y":[-0.5],"name":"b3caa38189be2e93dbbaeec4345e5d5f","type":"scatter","mode":"markers","opacity":0,"hoverinfo":"skip","showlegend":false,"marker":{"color":[0,1],"colorscale":[[0,"#132B43"],[0.0033444816053511323,"#132B44"],[0.0066889632107022647,"#132C44"],[0.010033444816053561,"#142C45"],[0.013377926421404692,"#142D45"],[0.016722408026755824,"#142D46"],[0.020066889632106958,"#142D46"],[0.023411371237458255,"#142E47"],[0.026755852842809385,"#152E47"],[0.030100334448160519,"#152F48"],[0.033444816053511649,"#152F48"],[0.036789297658862949,"#152F49"],[0.040133779264214076,"#153049"],[0.043478260869565209,"#16304A"],[0.046822742474916343,"#16304A"],[0.050167224080267636,"#16314B"],[0.05351170568561877,"#16314B"],[0.056856187290969903,"#16324C"],[0.060200668896321037,"#17324D"],[0.063545150501672171,"#17324D"],[0.066889632107023464,"#17334E"],[0.070234113712374591,"#17334E"],[0.073578595317725731,"#17344F"],[0.076923076923076858,"#18344F"],[0.080267558528428151,"#183450"],[0.083612040133779292,"#183550"],[0.086956521739130418,"#183551"],[0.090301003344481559,"#183651"],[0.093645484949832686,"#193652"],[0.096989966555183979,"#193652"],[0.10033444816053512,"#193753"],[0.10367892976588625,"#193754"],[0.10702341137123754,"#193854"],[0.11036789297658868,"#1A3855"],[0.11371237458193981,"#1A3955"],[0.11705685618729093,"#1A3956"],[0.12040133779264207,"#1A3956"],[0.12374581939799337,"#1A3A57"],[0.12709030100334451,"#1B3A57"],[0.13043478260869562,"#1B3B58"],[0.13377926421404676,"#1B3B59"],[0.13712374581939807,"#1B3B59"],[0.14046822742474918,"#1C3C5A"],[0.14381270903010032,"#1C3C5A"],[0.14715719063545146,"#1C3D5B"],[0.15050167224080258,"#1C3D5B"],[0.15384615384615388,"#1C3D5C"],[0.15719063545150502,"#1D3E5C"],[0.16053511705685614,"#1D3E5D"],[0.16387959866220744,"#1D3F5D"],[0.16722408026755858,"#1D3F5E"],[0.1705685618729097,"#1D3F5F"],[0.17391304347826084,"#1E405F"],[0.17725752508361198,"#1E4060"],[0.18060200668896326,"#1E4160"],[0.1839464882943144,"#1E4161"],[0.18729096989966554,"#1E4261"],[0.19063545150501668,"#1F4262"],[0.19397993311036796,"#1F4263"],[0.1973244147157191,"#1F4363"],[0.20066889632107024,"#1F4364"],[0.20401337792642135,"#1F4464"],[0.20735785953177249,"#204465"],[0.2107023411371238,"#204465"],[0.21404682274247491,"#204566"],[0.21739130434782605,"#204566"],[0.22073578595317736,"#214667"],[0.22408026755852847,"#214668"],[0.22742474916387961,"#214768"],[0.23076923076923075,"#214769"],[0.23411371237458187,"#214769"],[0.23745819397993317,"#22486A"],[0.24080267558528431,"#22486A"],[0.24414715719063543,"#22496B"],[0.24749163879598657,"#22496C"],[0.25083612040133785,"#224A6C"],[0.25418060200668902,"#234A6D"],[0.25752508361204013,"#234A6D"],[0.26086956521739124,"#234B6E"],[0.26421404682274241,"#234B6E"],[0.26755852842809369,"#244C6F"],[0.2709030100334448,"#244C70"],[0.27424749163879597,"#244C70"],[0.27759197324414725,"#244D71"],[0.28093645484949836,"#244D71"],[0.28428093645484953,"#254E72"],[0.28762541806020064,"#254E72"],[0.29096989966555176,"#254F73"],[0.29431438127090293,"#254F74"],[0.2976588628762542,"#254F74"],[0.30100334448160532,"#265075"],[0.30434782608695649,"#265075"],[0.30769230769230776,"#265176"],[0.31103678929765888,"#265176"],[0.31438127090301005,"#275277"],[0.31772575250836116,"#275278"],[0.32107023411371227,"#275278"],[0.32441471571906361,"#275379"],[0.32775919732441472,"#275379"],[0.33110367892976583,"#28547A"],[0.33444816053511717,"#28547B"],[0.33779264214046828,"#28557B"],[0.34113712374581939,"#28557C"],[0.34448160535117056,"#28567C"],[0.34782608695652167,"#29567D"],[0.35117056856187279,"#29567D"],[0.35451505016722412,"#29577E"],[0.35785953177257523,"#29577F"],[0.36120401337792635,"#2A587F"],[0.36454849498327768,"#2A5880"],[0.3678929765886288,"#2A5980"],[0.37123745819397991,"#2A5981"],[0.37458193979933108,"#2A5982"],[0.37792642140468219,"#2B5A82"],[0.38127090301003336,"#2B5A83"],[0.38461538461538464,"#2B5B83"],[0.38795986622073575,"#2B5B84"],[0.39130434782608703,"#2C5C85"],[0.3946488294314382,"#2C5C85"],[0.39799331103678931,"#2C5D86"],[0.40133779264214048,"#2C5D86"],[0.40468227424749159,"#2C5D87"],[0.4080267558528427,"#2D5E87"],[0.41137123745819404,"#2D5E88"],[0.41471571906354515,"#2D5F89"],[0.41806020066889643,"#2D5F89"],[0.4214046822742476,"#2E608A"],[0.42474916387959871,"#2E608A"],[0.42809364548494983,"#2E618B"],[0.43143812709030099,"#2E618C"],[0.43478260869565211,"#2E618C"],[0.43812709030100322,"#2F628D"],[0.44147157190635455,"#2F628D"],[0.44481605351170567,"#2F638E"],[0.44816053511705695,"#2F638F"],[0.45150501672240811,"#30648F"],[0.45484949832775923,"#306490"],[0.45819397993311034,"#306590"],[0.46153846153846151,"#306591"],[0.46488294314381262,"#306592"],[0.4682274247491639,"#316692"],[0.47157190635451507,"#316693"],[0.47491638795986635,"#316793"],[0.47826086956521746,"#316794"],[0.48160535117056863,"#326895"],[0.48494983277591974,"#326895"],[0.48829431438127086,"#326996"],[0.49163879598662202,"#326996"],[0.49498327759197314,"#326997"],[0.49832775919732442,"#336A98"],[0.50167224080267558,"#336A98"],[0.50501672240802686,"#336B99"],[0.50836120401337803,"#336B99"],[0.51170568561872909,"#346C9A"],[0.51505016722408026,"#346C9B"],[0.51839464882943143,"#346D9B"],[0.52173913043478248,"#346D9C"],[0.52508361204013387,"#346E9D"],[0.52842809364548493,"#356E9D"],[0.53177257525083621,"#356E9E"],[0.53511705685618738,"#356F9E"],[0.53846153846153855,"#356F9F"],[0.5418060200668896,"#3670A0"],[0.54515050167224077,"#3670A0"],[0.54849498327759194,"#3671A1"],[0.551839464882943,"#3671A1"],[0.55518394648829417,"#3772A2"],[0.55852842809364533,"#3772A3"],[0.56187290969899673,"#3773A3"],[0.56521739130434789,"#3773A4"],[0.56856187290969906,"#3773A4"],[0.57190635451505012,"#3874A5"],[0.57525083612040129,"#3874A6"],[0.57859531772575246,"#3875A6"],[0.58193979933110385,"#3875A7"],[0.58528428093645501,"#3976A8"],[0.58862876254180618,"#3976A8"],[0.59197324414715724,"#3977A9"],[0.59531772575250841,"#3977A9"],[0.59866220735785958,"#3978AA"],[0.60200668896321063,"#3A78AB"],[0.6053511705685618,"#3A79AB"],[0.60869565217391297,"#3A79AC"],[0.61204013377926403,"#3A79AC"],[0.6153846153846152,"#3B7AAD"],[0.6187290969899667,"#3B7AAE"],[0.62207357859531776,"#3B7BAE"],[0.62541806020066892,"#3B7BAF"],[0.62876254180602009,"#3C7CB0"],[0.63210702341137115,"#3C7CB0"],[0.63545150501672232,"#3C7DB1"],[0.63879598662207382,"#3C7DB1"],[0.64214046822742488,"#3C7EB2"],[0.64548494983277604,"#3D7EB3"],[0.64882943143812721,"#3D7FB3"],[0.65217391304347827,"#3D7FB4"],[0.65551839464882944,"#3D7FB5"],[0.65886287625418061,"#3E80B5"],[0.66220735785953166,"#3E80B6"],[0.66555183946488283,"#3E81B6"],[0.668896321070234,"#3E81B7"],[0.67224080267558506,"#3F82B8"],[0.67558528428093656,"#3F82B8"],[0.67892976588628773,"#3F83B9"],[0.68227424749163879,"#3F83BA"],[0.68561872909698995,"#4084BA"],[0.68896321070234112,"#4084BB"],[0.69230769230769218,"#4085BB"],[0.69565217391304368,"#4085BC"],[0.69899665551839485,"#4086BD"],[0.70234113712374591,"#4186BD"],[0.70568561872909707,"#4186BE"],[0.70903010033444824,"#4187BF"],[0.7123745819397993,"#4187BF"],[0.71571906354515047,"#4288C0"],[0.71906354515050164,"#4288C1"],[0.72240802675585269,"#4289C1"],[0.72575250836120386,"#4289C2"],[0.72909698996655503,"#438AC2"],[0.73244147157190642,"#438AC3"],[0.73578595317725759,"#438BC4"],[0.73913043478260876,"#438BC4"],[0.74247491638795982,"#438CC5"],[0.74581939799331098,"#448CC6"],[0.74916387959866215,"#448DC6"],[0.75250836120401354,"#448DC7"],[0.75585284280936471,"#448EC8"],[0.75919732441471588,"#458EC8"],[0.76254180602006694,"#458FC9"],[0.76588628762541811,"#458FC9"],[0.76923076923076927,"#458FCA"],[0.77257525083612044,"#4690CB"],[0.7759197324414715,"#4690CB"],[0.77926421404682267,"#4691CC"],[0.78260869565217384,"#4691CD"],[0.78595317725752489,"#4792CD"],[0.78929765886287606,"#4792CE"],[0.79264214046822756,"#4793CF"],[0.79598662207357862,"#4793CF"],[0.79933110367892979,"#4894D0"],[0.80267558528428096,"#4894D0"],[0.80602006688963201,"#4895D1"],[0.80936454849498352,"#4895D2"],[0.81270903010033468,"#4896D2"],[0.81605351170568574,"#4996D3"],[0.81939799331103691,"#4997D4"],[0.82274247491638808,"#4997D4"],[0.82608695652173914,"#4998D5"],[0.8294314381270903,"#4A98D6"],[0.83277591973244147,"#4A99D6"],[0.83612040133779253,"#4A99D7"],[0.8394648829431437,"#4A9AD8"],[0.84280936454849487,"#4B9AD8"],[0.84615384615384592,"#4B9BD9"],[0.84949832775919742,"#4B9BDA"],[0.85284280936454859,"#4B9BDA"],[0.85618729096989965,"#4C9CDB"],[0.85953177257525082,"#4C9CDB"],[0.86287625418060199,"#4C9DDC"],[0.86622073578595338,"#4C9DDD"],[0.86956521739130455,"#4D9EDD"],[0.87290969899665571,"#4D9EDE"],[0.87625418060200677,"#4D9FDF"],[0.87959866220735794,"#4D9FDF"],[0.88294314381270911,"#4DA0E0"],[0.88628762541806017,"#4EA0E1"],[0.88963210702341133,"#4EA1E1"],[0.8929765886287625,"#4EA1E2"],[0.89632107023411356,"#4EA2E3"],[0.89966555183946473,"#4FA2E3"],[0.9030100334448159,"#4FA3E4"],[0.90635451505016729,"#4FA3E5"],[0.90969899665551845,"#4FA4E5"],[0.91304347826086962,"#50A4E6"],[0.91638795986622068,"#50A5E7"],[0.91973244147157185,"#50A5E7"],[0.92307692307692335,"#50A6E8"],[0.92642140468227441,"#51A6E8"],[0.92976588628762558,"#51A7E9"],[0.93311036789297674,"#51A7EA"],[0.9364548494983278,"#51A8EA"],[0.93979933110367897,"#52A8EB"],[0.94314381270903014,"#52A9EC"],[0.9464882943143812,"#52A9EC"],[0.94983277591973236,"#52AAED"],[0.95317725752508353,"#53AAEE"],[0.95652173913043459,"#53ABEE"],[0.95986622073578576,"#53ABEF"],[0.96321070234113726,"#53ACF0"],[0.96655518394648832,"#54ACF0"],[0.96989966555183948,"#54ADF1"],[0.97324414715719065,"#54ADF2"],[0.97658862876254171,"#54AEF2"],[0.97993311036789321,"#55AEF3"],[0.98327759197324438,"#55AFF4"],[0.98662207357859544,"#55AFF4"],[0.98996655518394661,"#55B0F5"],[0.99331103678929777,"#56B0F6"],[0.99665551839464883,"#56B1F6"],[1,"#56B1F7"]],"colorbar":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"thickness":23.040000000000003,"title":"colData(se.subset)$TIN_median","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"tickmode":"array","ticktext":["50","60","70","80"],"tickvals":[0.22724018475750571,0.45741724403387213,0.68759430331023863,0.91777136258660508],"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"ticklen":2,"len":0.5,"yanchor":"top","y":1}},"xaxis":"x","yaxis":"y","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":37.260273972602747,"l":54.794520547945211},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.45062006448861935,0.41842937640520583],"tickmode":"array","ticktext":["-0.4","-0.2","0.0","0.2","0.4"],"tickvals":[-0.40000000000000002,-0.20000000000000001,0,0.20000000000000007,0.40000000000000002],"categoryorder":"array","categoryarray":["-0.4","-0.2","0.0","0.2","0.4"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":{"text":"PC1: 28% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.6665199896886409,0.46067339298057014],"tickmode":"array","ticktext":["-0.50","-0.25","0.00","0.25"],"tickvals":[-0.5,-0.25,0,0.25],"categoryorder":"array","categoryarray":["-0.50","-0.25","0.00","0.25"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"PC2: 18% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"y":0.5,"yanchor":"top","title":{"text":"colData(se.subset)$collection_site<br />colData(se.subset)$TIN_median","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"692055e68d52":{"x":{},"y":{},"colour":{},"label":{},"shape":{},"type":"scatter"}},"cur_data":"692055e68d52","visdat":{"692055e68d52":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```


The PCA plot shows the first two principal components of the log cpm transformed expression data 
TIN values illustrates for the variation that is observed and shape indicates where samples originate. 

### PCA plot colored by RIN 


```r
#ggplot(as.data.frame(pca$rotation[,c(1,2)]),aes(PC1,PC2, color=colData(se)$RIN, shape=colData(se)$tissue_type)) +
plot3 <- ggplot(as.data.frame(pca$rotation[,c(1,2)]),
                aes(PC1,PC2, color=as.numeric(colData(se.subset)$RIN),label = sample_id)) +
  geom_point(size=3, aes(shape=colData(se.subset)$collection_site)) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
#plot3
ggplotly(plot3)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-384a0bf4dc3dd7668aec" style="width:1056px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-384a0bf4dc3dd7668aec">{"x":{"data":[{"x":[-0.41111781717526369,0.3406716837676157,-0.22213877826446352,0.24694489205281545,-0.39487515083993557,0.32717429264698528,-0.23990188662012524,0.25188371136087773,-0.27756807602035954,0.37892712909185017],"y":[0.40943733013196965,0.11972671630713563,0.20647913029965592,0.19226077096856722,-0.28166110669463779,-0.049867851015437524,-0.61528392684004041,-0.42020469336764826,0.26981755589983103,0.16929607431060212],"text":["PC1: -0.4111178<br />PC2:  0.40943733<br />as.numeric(colData(se.subset)$RIN): 9.0<br />sample_id: SynCohort_P01_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.3406717<br />PC2:  0.11972672<br />as.numeric(colData(se.subset)$RIN): 9.2<br />sample_id: SynCohort_P02_S1<br />colData(se.subset)$collection_site: SiteA","PC1: -0.2221388<br />PC2:  0.20647913<br />as.numeric(colData(se.subset)$RIN): 5.6<br />sample_id: SynCohort_P03_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.2469449<br />PC2:  0.19226077<br />as.numeric(colData(se.subset)$RIN): 8.6<br />sample_id: SynCohort_P04_S1<br />colData(se.subset)$collection_site: SiteA","PC1: -0.3948752<br />PC2: -0.28166111<br />as.numeric(colData(se.subset)$RIN): 7.5<br />sample_id: SynCohort_P05_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.3271743<br />PC2: -0.04986785<br />as.numeric(colData(se.subset)$RIN): 6.9<br />sample_id: SynCohort_P06_S1<br />colData(se.subset)$collection_site: SiteA","PC1: -0.2399019<br />PC2: -0.61528393<br />as.numeric(colData(se.subset)$RIN): 8.1<br />sample_id: SynCohort_P07_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.2518837<br />PC2: -0.42020469<br />as.numeric(colData(se.subset)$RIN): 4.7<br />sample_id: SynCohort_P08_S1<br />colData(se.subset)$collection_site: SiteA","PC1: -0.2775681<br />PC2:  0.26981756<br />as.numeric(colData(se.subset)$RIN): 7.6<br />sample_id: SynCohort_P09_S1<br />colData(se.subset)$collection_site: SiteA","PC1:  0.3789271<br />PC2:  0.16929607<br />as.numeric(colData(se.subset)$RIN): 7.9<br />sample_id: SynCohort_P10_S1<br />colData(se.subset)$collection_site: SiteA"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(83,170,238,1)","rgba(86,177,247,1)","rgba(31,67,100,1)","rgba(76,157,221,1)","rgba(59,123,174,1)","rgba(50,105,150,1)","rgba(68,142,199,1)","rgba(19,43,67,1)","rgba(61,126,179,1)","rgba(65,135,191,1)"],"opacity":1,"size":11.338582677165356,"symbol":"circle","line":{"width":1.8897637795275593,"color":["rgba(83,170,238,1)","rgba(86,177,247,1)","rgba(31,67,100,1)","rgba(76,157,221,1)","rgba(59,123,174,1)","rgba(50,105,150,1)","rgba(68,142,199,1)","rgba(19,43,67,1)","rgba(61,126,179,1)","rgba(65,135,191,1)"]}},"hoveron":"points","name":"SiteA","legendgroup":"SiteA","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.40000000000000002],"y":[-0.5],"name":"ec566e40515098847ca6946ab34ef1db","type":"scatter","mode":"markers","opacity":0,"hoverinfo":"skip","showlegend":false,"marker":{"color":[0,1],"colorscale":[[0,"#132B43"],[0.0033444816053511935,"#132B44"],[0.006688963210702387,"#132C44"],[0.010033444816053581,"#142C45"],[0.013377926421404774,"#142D45"],[0.016722408026755772,"#142D46"],[0.020066889632106965,"#142D46"],[0.023411371237458158,"#142E47"],[0.026755852842809354,"#152E47"],[0.030100334448160546,"#152F48"],[0.033444816053511739,"#152F48"],[0.036789297658862935,"#152F49"],[0.040133779264214124,"#153049"],[0.043478260869565126,"#16304A"],[0.046822742474916315,"#16304A"],[0.050167224080267511,"#16314B"],[0.053511705685618707,"#16314B"],[0.056856187290969896,"#16324C"],[0.060200668896321093,"#17324D"],[0.063545150501672282,"#17324D"],[0.066889632107023478,"#17334E"],[0.070234113712374674,"#17334E"],[0.073578595317725662,"#17344F"],[0.076923076923076858,"#18344F"],[0.080267558528428054,"#183450"],[0.08361204013377925,"#183550"],[0.086956521739130446,"#183551"],[0.090301003344481642,"#183651"],[0.093645484949832825,"#193652"],[0.096989966555184021,"#193652"],[0.10033444816053522,"#193753"],[0.10367892976588622,"#193754"],[0.10702341137123741,"#193854"],[0.1103678929765886,"#1A3855"],[0.11371237458193979,"#1A3955"],[0.11705685618729099,"#1A3956"],[0.12040133779264219,"#1A3956"],[0.12374581939799338,"#1A3A57"],[0.12709030100334437,"#1B3A57"],[0.13043478260869557,"#1B3B58"],[0.13377926421404676,"#1B3B59"],[0.13712374581939796,"#1B3B59"],[0.14046822742474915,"#1C3C5A"],[0.14381270903010035,"#1C3C5A"],[0.14715719063545155,"#1C3D5B"],[0.15050167224080274,"#1C3D5B"],[0.15384615384615391,"#1C3D5C"],[0.15719063545150511,"#1D3E5C"],[0.16053511705685611,"#1D3E5D"],[0.1638795986622073,"#1D3F5D"],[0.1672240802675585,"#1D3F5E"],[0.1705685618729097,"#1D3F5F"],[0.17391304347826089,"#1E405F"],[0.17725752508361209,"#1E4060"],[0.18060200668896328,"#1E4160"],[0.18394648829431426,"#1E4161"],[0.18729096989966545,"#1E4261"],[0.19063545150501665,"#1F4262"],[0.19397993311036785,"#1F4263"],[0.19732441471571904,"#1F4363"],[0.20066889632107024,"#1F4364"],[0.20401337792642144,"#1F4464"],[0.20735785953177263,"#204465"],[0.21070234113712383,"#204465"],[0.21404682274247502,"#204566"],[0.217391304347826,"#204566"],[0.22073578595317719,"#214667"],[0.22408026755852839,"#214668"],[0.22742474916387959,"#214768"],[0.23076923076923078,"#214769"],[0.23411371237458198,"#214769"],[0.23745819397993317,"#22486A"],[0.24080267558528418,"#22486A"],[0.24414715719063537,"#22496B"],[0.24749163879598657,"#22496C"],[0.25083612040133774,"#224A6C"],[0.25418060200668896,"#234A6D"],[0.25752508361204013,"#234A6D"],[0.26086956521739135,"#234B6E"],[0.26421404682274252,"#234B6E"],[0.26755852842809374,"#244C6F"],[0.27090301003344491,"#244C70"],[0.27424749163879608,"#244C70"],[0.27759197324414708,"#244D71"],[0.28093645484949831,"#244D71"],[0.28428093645484948,"#254E72"],[0.2876254180602007,"#254E72"],[0.29096989966555187,"#254F73"],[0.29431438127090309,"#254F74"],[0.29765886287625409,"#254F74"],[0.30100334448160526,"#265075"],[0.30434782608695643,"#265075"],[0.30769230769230765,"#265176"],[0.31103678929765882,"#265176"],[0.31438127090301005,"#275277"],[0.31772575250836121,"#275278"],[0.32107023411371244,"#275278"],[0.32441471571906361,"#275379"],[0.32775919732441483,"#275379"],[0.331103678929766,"#28547A"],[0.334448160535117,"#28547B"],[0.33779264214046817,"#28557B"],[0.34113712374581939,"#28557C"],[0.34448160535117056,"#28567C"],[0.34782608695652178,"#29567D"],[0.35117056856187295,"#29567D"],[0.35451505016722396,"#29577E"],[0.35785953177257518,"#29577F"],[0.36120401337792635,"#2A587F"],[0.36454849498327757,"#2A5880"],[0.36789297658862874,"#2A5980"],[0.37123745819397991,"#2A5981"],[0.37458193979933113,"#2A5982"],[0.3779264214046823,"#2B5A82"],[0.38127090301003352,"#2B5A83"],[0.38461538461538469,"#2B5B83"],[0.38795986622073592,"#2B5B84"],[0.39130434782608692,"#2C5C85"],[0.39464882943143809,"#2C5C85"],[0.39799331103678931,"#2C5D86"],[0.40133779264214048,"#2C5D86"],[0.40468227424749165,"#2C5D87"],[0.40802675585284287,"#2D5E87"],[0.41137123745819387,"#2D5E88"],[0.41471571906354504,"#2D5F89"],[0.41806020066889626,"#2D5F89"],[0.42140468227424743,"#2E608A"],[0.42474916387959866,"#2E608A"],[0.42809364548494983,"#2E618B"],[0.43143812709030105,"#2E618C"],[0.43478260869565222,"#2E618C"],[0.43812709030100339,"#2F628D"],[0.44147157190635461,"#2F628D"],[0.44481605351170578,"#2F638E"],[0.448160535117057,"#2F638F"],[0.451505016722408,"#30648F"],[0.45484949832775917,"#306490"],[0.4581939799331104,"#306590"],[0.46153846153846156,"#306591"],[0.46488294314381257,"#306592"],[0.46822742474916373,"#316692"],[0.47157190635451496,"#316693"],[0.47491638795986613,"#316793"],[0.47826086956521735,"#316794"],[0.48160535117056852,"#326895"],[0.48494983277591974,"#326895"],[0.48829431438127091,"#326996"],[0.49163879598662213,"#326996"],[0.4949832775919733,"#326997"],[0.49832775919732453,"#336A98"],[0.5016722408026757,"#336A98"],[0.50501672240802686,"#336B99"],[0.50836120401337792,"#336B99"],[0.51170568561872909,"#346C9A"],[0.51505016722408026,"#346C9B"],[0.51839464882943143,"#346D9B"],[0.52173913043478248,"#346D9C"],[0.52508361204013365,"#346E9D"],[0.52842809364548482,"#356E9D"],[0.5317725752508361,"#356E9E"],[0.53511705685618727,"#356F9E"],[0.53846153846153844,"#356F9F"],[0.5418060200668896,"#3670A0"],[0.54515050167224077,"#3670A0"],[0.54849498327759205,"#3671A1"],[0.55183946488294322,"#3671A1"],[0.55518394648829439,"#3772A2"],[0.55852842809364556,"#3772A3"],[0.56187290969899684,"#3773A3"],[0.56521739130434778,"#3773A4"],[0.56856187290969895,"#3773A4"],[0.57190635451505012,"#3874A5"],[0.5752508361204014,"#3874A6"],[0.57859531772575234,"#3875A6"],[0.58193979933110351,"#3875A7"],[0.58528428093645479,"#3976A8"],[0.58862876254180596,"#3976A8"],[0.59197324414715713,"#3977A9"],[0.5953177257525083,"#3977A9"],[0.59866220735785958,"#3978AA"],[0.60200668896321075,"#3A78AB"],[0.60535117056856191,"#3A79AB"],[0.60869565217391308,"#3A79AC"],[0.61204013377926425,"#3A79AC"],[0.61538461538461553,"#3B7AAD"],[0.6187290969899667,"#3B7AAE"],[0.62207357859531764,"#3B7BAE"],[0.62541806020066892,"#3B7BAF"],[0.62876254180602009,"#3C7CB0"],[0.63210702341137126,"#3C7CB0"],[0.63545150501672221,"#3C7DB1"],[0.63879598662207349,"#3C7DB1"],[0.64214046822742465,"#3C7EB2"],[0.64548494983277582,"#3D7EB3"],[0.64882943143812699,"#3D7FB3"],[0.65217391304347827,"#3D7FB4"],[0.65551839464882944,"#3D7FB5"],[0.65886287625418061,"#3E80B5"],[0.66220735785953178,"#3E80B6"],[0.66555183946488305,"#3E81B6"],[0.66889632107023422,"#3E81B7"],[0.67224080267558539,"#3F82B8"],[0.67558528428093656,"#3F82B8"],[0.67892976588628762,"#3F83B9"],[0.68227424749163879,"#3F83BA"],[0.68561872909698995,"#4084BA"],[0.68896321070234112,"#4084BB"],[0.69230769230769218,"#4085BB"],[0.69565217391304335,"#4085BC"],[0.69899665551839452,"#4086BD"],[0.70234113712374568,"#4186BD"],[0.70568561872909696,"#4186BE"],[0.70903010033444813,"#4187BF"],[0.7123745819397993,"#4187BF"],[0.71571906354515047,"#4288C0"],[0.71906354515050175,"#4288C1"],[0.72240802675585292,"#4289C1"],[0.72575250836120409,"#4289C2"],[0.72909698996655525,"#438AC2"],[0.73244147157190653,"#438AC3"],[0.7357859531772577,"#438BC4"],[0.73913043478260843,"#438BC4"],[0.74247491638796004,"#438CC5"],[0.74581939799331087,"#448CC6"],[0.74916387959866249,"#448DC6"],[0.75250836120401321,"#448DC7"],[0.75585284280936449,"#448EC8"],[0.75919732441471566,"#458EC8"],[0.76254180602006683,"#458FC9"],[0.76588628762541799,"#458FC9"],[0.76923076923076916,"#458FCA"],[0.77257525083612044,"#4690CB"],[0.77591973244147161,"#4690CB"],[0.77926421404682278,"#4691CC"],[0.78260869565217395,"#4691CD"],[0.78595317725752523,"#4792CD"],[0.78929765886287639,"#4792CE"],[0.79264214046822756,"#4793CF"],[0.7959866220735784,"#4793CF"],[0.79933110367893001,"#4894D0"],[0.80267558528428073,"#4894D0"],[0.80602006688963235,"#4895D1"],[0.80936454849498318,"#4895D2"],[0.81270903010033435,"#4896D2"],[0.81605351170568552,"#4996D3"],[0.81939799331103669,"#4997D4"],[0.82274247491638797,"#4997D4"],[0.82608695652173914,"#4998D5"],[0.8294314381270903,"#4A98D6"],[0.83277591973244147,"#4A99D6"],[0.83612040133779264,"#4A99D7"],[0.83946488294314392,"#4A9AD8"],[0.84280936454849509,"#4B9AD8"],[0.84615384615384626,"#4B9BD9"],[0.84949832775919742,"#4B9BDA"],[0.85284280936454826,"#4B9BDA"],[0.85618729096989987,"#4C9CDB"],[0.85953177257525071,"#4C9CDB"],[0.86287625418060221,"#4C9DDC"],[0.86622073578595304,"#4C9DDD"],[0.86956521739130421,"#4D9EDD"],[0.87290969899665538,"#4D9EDE"],[0.87625418060200666,"#4D9FDF"],[0.87959866220735783,"#4D9FDF"],[0.882943143812709,"#4DA0E0"],[0.88628762541806017,"#4EA0E1"],[0.88963210702341144,"#4EA1E1"],[0.89297658862876261,"#4EA1E2"],[0.89632107023411378,"#4EA2E3"],[0.89966555183946495,"#4FA2E3"],[0.90301003344481579,"#4FA3E4"],[0.9063545150501674,"#4FA3E5"],[0.90969899665551812,"#4FA4E5"],[0.91304347826086973,"#50A4E6"],[0.91638795986622057,"#50A5E7"],[0.91973244147157218,"#50A5E7"],[0.92307692307692291,"#50A6E8"],[0.92642140468227452,"#51A6E8"],[0.92976588628762535,"#51A7E9"],[0.93311036789297652,"#51A7EA"],[0.93645484949832769,"#51A8EA"],[0.93979933110367886,"#52A8EB"],[0.94314381270903014,"#52A9EC"],[0.94648829431438131,"#52A9EC"],[0.94983277591973247,"#52AAED"],[0.95317725752508364,"#53AAEE"],[0.95652173913043492,"#53ABEE"],[0.95986622073578565,"#53ABEF"],[0.96321070234113726,"#53ACF0"],[0.96655518394648809,"#54ACF0"],[0.96989966555183971,"#54ADF1"],[0.97324414715719043,"#54ADF2"],[0.97658862876254204,"#54AEF2"],[0.97993311036789288,"#55AEF3"],[0.98327759197324438,"#55AFF4"],[0.98662207357859522,"#55AFF4"],[0.98996655518394638,"#55B0F5"],[0.99331103678929766,"#56B0F6"],[0.99665551839464883,"#56B1F6"],[1,"#56B1F7"]],"colorbar":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"thickness":23.040000000000003,"title":"as.numeric(colData(se.subset)$RIN)","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"tickmode":"array","ticktext":["5","6","7","8","9"],"tickvals":[0.068111111111111081,0.28959259259259257,0.51107407407407413,0.73255555555555563,0.95403703703703713],"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"ticklen":2,"len":0.5,"yanchor":"top","y":1}},"xaxis":"x","yaxis":"y","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":37.260273972602747,"l":54.794520547945211},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.45062006448861935,0.41842937640520583],"tickmode":"array","ticktext":["-0.4","-0.2","0.0","0.2","0.4"],"tickvals":[-0.40000000000000002,-0.20000000000000001,0,0.20000000000000007,0.40000000000000002],"categoryorder":"array","categoryarray":["-0.4","-0.2","0.0","0.2","0.4"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"y","title":{"text":"PC1: 28% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.6665199896886409,0.46067339298057014],"tickmode":"array","ticktext":["-0.50","-0.25","0.00","0.25"],"tickvals":[-0.5,-0.25,0,0.25],"categoryorder":"array","categoryarray":["-0.50","-0.25","0.00","0.25"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.68949771689498},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176002,"zeroline":false,"anchor":"x","title":{"text":"PC2: 18% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}},"hoverformat":".2f"},"shapes":[],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.68949771689498},"y":0.5,"yanchor":"top","title":{"text":"colData(se.subset)$collection_site<br />as.numeric(colData(se.subset)$RIN)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"6920cd2ffab":{"x":{},"y":{},"colour":{},"label":{},"shape":{},"type":"scatter"}},"cur_data":"6920cd2ffab","visdat":{"6920cd2ffab":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```


The PCA plot shows the first two principal components of the log cpm transformed expression data 
RIN values illustrates for the variation that is observed and shape indicates where samples originate. 

### PCA plot colored by outliers 


```r
#ggplot(as.data.frame(pca$rotation[,c(1,2)]), aes(PC1,PC2, color=out.color, shape=colData(se)$tissue_type)) +
plot4 <- ggplot(as.data.frame(pca$rotation[,c(1,2)]), 
                aes(PC1,PC2, shape=colData(se.subset)$collection_site, label = sample_id)) +
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
  labs(color='Outliers',shape= "Site")+
  ggplot2::scale_color_manual(name = "Outliers",
                     values = c("firebrick", "blue"),
                     labels = c("no", "yes"))
  #labs(color='Outliers',shape= "Site") 
  

ggplotly(plot4)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-03fd4504c77dab3f9e87" style="width:960px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-03fd4504c77dab3f9e87">{"x":{"data":[{"x":[-0.41111781717526369,0.3406716837676157,-0.22213877826446352,0.24694489205281545,-0.39487515083993557,0.32717429264698528,-0.23990188662012524,0.25188371136087773,-0.27756807602035954,0.37892712909185017],"y":[0.40943733013196965,0.11972671630713563,0.20647913029965592,0.19226077096856722,-0.28166110669463779,-0.049867851015437524,-0.61528392684004041,-0.42020469336764826,0.26981755589983103,0.16929607431060212],"text":["PC1: -0.4111178<br />PC2:  0.40943733<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P01_S1<br />out.color: no","PC1:  0.3406717<br />PC2:  0.11972672<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P02_S1<br />out.color: no","PC1: -0.2221388<br />PC2:  0.20647913<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P03_S1<br />out.color: no","PC1:  0.2469449<br />PC2:  0.19226077<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P04_S1<br />out.color: no","PC1: -0.3948752<br />PC2: -0.28166111<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P05_S1<br />out.color: no","PC1:  0.3271743<br />PC2: -0.04986785<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P06_S1<br />out.color: no","PC1: -0.2399019<br />PC2: -0.61528393<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P07_S1<br />out.color: no","PC1:  0.2518837<br />PC2: -0.42020469<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P08_S1<br />out.color: no","PC1: -0.2775681<br />PC2:  0.26981756<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P09_S1<br />out.color: no","PC1:  0.3789271<br />PC2:  0.16929607<br />colData(se.subset)$collection_site: SiteA<br />sample_id: SynCohort_P10_S1<br />out.color: no"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(178,34,34,1)","opacity":1,"size":15.118110236220474,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(178,34,34,1)"}},"hoveron":"points","name":"(no,SiteA)","legendgroup":"(no,SiteA)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":48.152760481527615,"l":76.048152760481514},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.45062006448861935,0.41842937640520583],"tickmode":"array","ticktext":["-0.4","-0.2","0.0","0.2","0.4"],"tickvals":[-0.40000000000000002,-0.20000000000000001,0,0.20000000000000007,0.40000000000000002],"categoryorder":"array","categoryarray":["-0.4","-0.2","0.0","0.2","0.4"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":18.596928185969286},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"PC1 (28%)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.6665199896886409,0.46067339298057014],"tickmode":"array","ticktext":["-0.50","-0.25","0.00","0.25"],"tickvals":[-0.5,-0.25,0,0.25],"categoryorder":"array","categoryarray":["-0.50","-0.25","0.00","0.25"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":18.596928185969279},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"PC2 (18%)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"rgba(255,255,255,1)","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":17.268576172685766},"title":{"text":"Site<br />Outliers","font":{"color":"rgba(0,0,0,1)","family":"","size":19.925280199252807}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"6920662b57d":{"x":{},"y":{},"shape":{},"label":{},"colour":{},"type":"scatter"}},"cur_data":"6920662b57d","visdat":{"6920662b57d":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
batch.plot <- ggplot(as.data.frame(pca.before$rotation[,c(1,2)]), 
                     aes(PC1,PC2, color=batch.sub, shape=colData(se.subset)$tissue_type,label=sampleID)) +
  geom_point(size=3) +
  labs(#title = "Before batch effect correction",
       x = paste0("PC1: ",percentVar[1],"% variance"),
       y = paste0("PC2: ",percentVar[2],"% variance")) +
   scale_shape_manual(name = "Sample Type",
                       values = 16,
                      labels =  tissue) +
  scale_color_manual(name = tissue,
                     #values = COLORS
                     values = myCol)+
  theme_classic() +
  theme(text = element_text(size = 14)) 


ggplotly(batch.plot)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-3b8232589d3760fce280" style="width:960px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-3b8232589d3760fce280">{"x":{"data":[{"x":[-0.41111781717526369,-0.22213877826446352,-0.39487515083993557,-0.23990188662012524,-0.27756807602035954],"y":[0.40943733013196965,0.20647913029965592,-0.28166110669463779,-0.61528392684004041,0.26981755589983103],"text":["PC1: -0.4111178<br />PC2:  0.40943733<br />batch.sub: Batch1<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P01_S1","PC1: -0.2221388<br />PC2:  0.20647913<br />batch.sub: Batch1<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P03_S1","PC1: -0.3948752<br />PC2: -0.28166111<br />batch.sub: Batch1<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P05_S1","PC1: -0.2399019<br />PC2: -0.61528393<br />batch.sub: Batch1<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P07_S1","PC1: -0.2775681<br />PC2:  0.26981756<br />batch.sub: Batch1<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P09_S1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(166,206,227,1)","opacity":1,"size":11.338582677165356,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(166,206,227,1)"}},"hoveron":"points","name":"(Lung,Batch1)","legendgroup":"(Lung,Batch1)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[0.3406716837676157,0.24694489205281545,0.32717429264698528,0.25188371136087773,0.37892712909185017],"y":[0.11972671630713563,0.19226077096856722,-0.049867851015437524,-0.42020469336764826,0.16929607431060212],"text":["PC1:  0.3406717<br />PC2:  0.11972672<br />batch.sub: Batch2<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P02_S1","PC1:  0.2469449<br />PC2:  0.19226077<br />batch.sub: Batch2<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P04_S1","PC1:  0.3271743<br />PC2: -0.04986785<br />batch.sub: Batch2<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P06_S1","PC1:  0.2518837<br />PC2: -0.42020469<br />batch.sub: Batch2<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P08_S1","PC1:  0.3789271<br />PC2:  0.16929607<br />batch.sub: Batch2<br />colData(se.subset)$tissue_type: Lung<br />sampleID: SynCohort_P10_S1"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":"rgba(177,89,40,1)","opacity":1,"size":11.338582677165356,"symbol":"circle","line":{"width":1.8897637795275593,"color":"rgba(177,89,40,1)"}},"hoveron":"points","name":"(Lung,Batch2)","legendgroup":"(Lung,Batch2)","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":44.433374844333756,"l":66.749688667496883},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.45062006448861935,0.41842937640520583],"tickmode":"array","ticktext":["-0.4","-0.2","0.0","0.2","0.4"],"tickvals":[-0.40000000000000002,-0.20000000000000001,0,0.20000000000000007,0.40000000000000002],"categoryorder":"array","categoryarray":["-0.4","-0.2","0.0","0.2","0.4"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":14.87754254877543},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176002,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"PC1: 28% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.6665199896886409,0.46067339298057014],"tickmode":"array","ticktext":["-0.50","-0.25","0.00","0.25"],"tickvals":[-0.5,-0.25,0,0.25],"categoryorder":"array","categoryarray":["-0.50","-0.25","0.00","0.25"],"nticks":null,"ticks":"outside","tickcolor":"rgba(0,0,0,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":14.877542548775427},"tickangle":-0,"showline":true,"linecolor":"rgba(0,0,0,1)","linewidth":0.66417600664176002,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"PC2: 18% variance","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"shapes":[],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":14.87754254877543},"title":{"text":"Sample Type<br />Lung","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"6920530264ee":{"x":{},"y":{},"colour":{},"shape":{},"label":{},"type":"scatter"}},"cur_data":"6920530264ee","visdat":{"6920530264ee":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
# NULL, i.e. when no QC metric happens to be constant across this tissue/site subset
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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/histo_tally_outliers-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/uniquely_mapped_percent-1.png" style="display: block; margin: auto;" />

Correlation between RSeQC unique mapped percent and STAR uniquely mapped pecent  colored by number of outliers binned.


### Correlaton of RIN vs TIN 


```r
ggscatter(qc.df2, x = "RIN", y = "TIN_median", color = "Number_Outliers_Binned",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "RNA integrity number (RIN)", ylab = "Transcript integrity number (TIN")
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/RIN_vs_TIN_median-1.png" style="display: block; margin: auto;" />

Correlation between Transcript integrity number (TIN) and RNA integrity number (RIN) colored by number of outliers binned.
 

### Correlation uniquelly mapped percent vs genes detected 


```r
qc.df2$Genes_detected <- colSums(assay(se.subset) > 0)
ggscatter(qc.df2, x = "STAR_uniquely_mapped_percent", y = "Genes_detected", color = "Number_Outliers_Binned",
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "STAR uniquely mapped percent", ylab = "Genes detected")
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/uniquely_mapped_percent_vs_gene_detected-1.png" style="display: block; margin: auto;" />

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
ggplotly(plotpca)
```

```{=html}
<div class="plotly html-widget html-fill-item" id="htmlwidget-7d0deb8db741c042f042" style="width:1152px;height:768px;"></div>
<script type="application/json" data-for="htmlwidget-7d0deb8db741c042f042">{"x":{"data":[{"x":[0.52780498846031543,-0.32593090929049179,0.29987422200967745,-0.43770760586155499,-0.25095646717192166,-0.037341457872466823,-0.28453662228253201,-0.081394999219723416,0.23362398148138813,0.3565648697473095],"y":[0.1018699799826475,0.47976097060731288,-0.33176018707094324,-0.34025063457475113,-0.48111972481703918,-0.25513948943024478,0.40931163800087522,0.20702332482216668,0.15277383956789886,0.057530282912078332],"text":["PC1:  0.52780499<br />PC2:  0.10186998<br />as.numeric(colData(se.subset)$RIN): 9.0<br />sampleID: SynCohort_P01_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.32593091<br />PC2:  0.47976097<br />as.numeric(colData(se.subset)$RIN): 9.2<br />sampleID: SynCohort_P02_S1<br />Number_Outliers_Binned: x =< 6","PC1:  0.29987422<br />PC2: -0.33176019<br />as.numeric(colData(se.subset)$RIN): 5.6<br />sampleID: SynCohort_P03_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.43770761<br />PC2: -0.34025063<br />as.numeric(colData(se.subset)$RIN): 8.6<br />sampleID: SynCohort_P04_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.25095647<br />PC2: -0.48111972<br />as.numeric(colData(se.subset)$RIN): 7.5<br />sampleID: SynCohort_P05_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.03734146<br />PC2: -0.25513949<br />as.numeric(colData(se.subset)$RIN): 6.9<br />sampleID: SynCohort_P06_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.28453662<br />PC2:  0.40931164<br />as.numeric(colData(se.subset)$RIN): 8.1<br />sampleID: SynCohort_P07_S1<br />Number_Outliers_Binned: x =< 6","PC1: -0.08139500<br />PC2:  0.20702332<br />as.numeric(colData(se.subset)$RIN): 4.7<br />sampleID: SynCohort_P08_S1<br />Number_Outliers_Binned: x =< 6","PC1:  0.23362398<br />PC2:  0.15277384<br />as.numeric(colData(se.subset)$RIN): 7.6<br />sampleID: SynCohort_P09_S1<br />Number_Outliers_Binned: x =< 6","PC1:  0.35656487<br />PC2:  0.05753028<br />as.numeric(colData(se.subset)$RIN): 7.9<br />sampleID: SynCohort_P10_S1<br />Number_Outliers_Binned: x =< 6"],"type":"scatter","mode":"markers","marker":{"autocolorscale":false,"color":["rgba(83,170,238,1)","rgba(86,177,247,1)","rgba(31,67,100,1)","rgba(76,157,221,1)","rgba(59,123,174,1)","rgba(50,105,150,1)","rgba(68,142,199,1)","rgba(19,43,67,1)","rgba(61,126,179,1)","rgba(65,135,191,1)"],"opacity":1,"size":18.897637795275593,"symbol":"circle","line":{"width":1.8897637795275593,"color":["rgba(83,170,238,1)","rgba(86,177,247,1)","rgba(31,67,100,1)","rgba(76,157,221,1)","rgba(59,123,174,1)","rgba(50,105,150,1)","rgba(68,142,199,1)","rgba(19,43,67,1)","rgba(61,126,179,1)","rgba(65,135,191,1)"]}},"hoveron":"points","name":"x =< 6","legendgroup":"x =< 6","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"x":[-0.25],"y":[-0.5],"name":"ec566e40515098847ca6946ab34ef1db","type":"scatter","mode":"markers","opacity":0,"hoverinfo":"skip","showlegend":false,"marker":{"color":[0,1],"colorscale":[[0,"#132B43"],[0.0033444816053511935,"#132B44"],[0.006688963210702387,"#132C44"],[0.010033444816053581,"#142C45"],[0.013377926421404774,"#142D45"],[0.016722408026755772,"#142D46"],[0.020066889632106965,"#142D46"],[0.023411371237458158,"#142E47"],[0.026755852842809354,"#152E47"],[0.030100334448160546,"#152F48"],[0.033444816053511739,"#152F48"],[0.036789297658862935,"#152F49"],[0.040133779264214124,"#153049"],[0.043478260869565126,"#16304A"],[0.046822742474916315,"#16304A"],[0.050167224080267511,"#16314B"],[0.053511705685618707,"#16314B"],[0.056856187290969896,"#16324C"],[0.060200668896321093,"#17324D"],[0.063545150501672282,"#17324D"],[0.066889632107023478,"#17334E"],[0.070234113712374674,"#17334E"],[0.073578595317725662,"#17344F"],[0.076923076923076858,"#18344F"],[0.080267558528428054,"#183450"],[0.08361204013377925,"#183550"],[0.086956521739130446,"#183551"],[0.090301003344481642,"#183651"],[0.093645484949832825,"#193652"],[0.096989966555184021,"#193652"],[0.10033444816053522,"#193753"],[0.10367892976588622,"#193754"],[0.10702341137123741,"#193854"],[0.1103678929765886,"#1A3855"],[0.11371237458193979,"#1A3955"],[0.11705685618729099,"#1A3956"],[0.12040133779264219,"#1A3956"],[0.12374581939799338,"#1A3A57"],[0.12709030100334437,"#1B3A57"],[0.13043478260869557,"#1B3B58"],[0.13377926421404676,"#1B3B59"],[0.13712374581939796,"#1B3B59"],[0.14046822742474915,"#1C3C5A"],[0.14381270903010035,"#1C3C5A"],[0.14715719063545155,"#1C3D5B"],[0.15050167224080274,"#1C3D5B"],[0.15384615384615391,"#1C3D5C"],[0.15719063545150511,"#1D3E5C"],[0.16053511705685611,"#1D3E5D"],[0.1638795986622073,"#1D3F5D"],[0.1672240802675585,"#1D3F5E"],[0.1705685618729097,"#1D3F5F"],[0.17391304347826089,"#1E405F"],[0.17725752508361209,"#1E4060"],[0.18060200668896328,"#1E4160"],[0.18394648829431426,"#1E4161"],[0.18729096989966545,"#1E4261"],[0.19063545150501665,"#1F4262"],[0.19397993311036785,"#1F4263"],[0.19732441471571904,"#1F4363"],[0.20066889632107024,"#1F4364"],[0.20401337792642144,"#1F4464"],[0.20735785953177263,"#204465"],[0.21070234113712383,"#204465"],[0.21404682274247502,"#204566"],[0.217391304347826,"#204566"],[0.22073578595317719,"#214667"],[0.22408026755852839,"#214668"],[0.22742474916387959,"#214768"],[0.23076923076923078,"#214769"],[0.23411371237458198,"#214769"],[0.23745819397993317,"#22486A"],[0.24080267558528418,"#22486A"],[0.24414715719063537,"#22496B"],[0.24749163879598657,"#22496C"],[0.25083612040133774,"#224A6C"],[0.25418060200668896,"#234A6D"],[0.25752508361204013,"#234A6D"],[0.26086956521739135,"#234B6E"],[0.26421404682274252,"#234B6E"],[0.26755852842809374,"#244C6F"],[0.27090301003344491,"#244C70"],[0.27424749163879608,"#244C70"],[0.27759197324414708,"#244D71"],[0.28093645484949831,"#244D71"],[0.28428093645484948,"#254E72"],[0.2876254180602007,"#254E72"],[0.29096989966555187,"#254F73"],[0.29431438127090309,"#254F74"],[0.29765886287625409,"#254F74"],[0.30100334448160526,"#265075"],[0.30434782608695643,"#265075"],[0.30769230769230765,"#265176"],[0.31103678929765882,"#265176"],[0.31438127090301005,"#275277"],[0.31772575250836121,"#275278"],[0.32107023411371244,"#275278"],[0.32441471571906361,"#275379"],[0.32775919732441483,"#275379"],[0.331103678929766,"#28547A"],[0.334448160535117,"#28547B"],[0.33779264214046817,"#28557B"],[0.34113712374581939,"#28557C"],[0.34448160535117056,"#28567C"],[0.34782608695652178,"#29567D"],[0.35117056856187295,"#29567D"],[0.35451505016722396,"#29577E"],[0.35785953177257518,"#29577F"],[0.36120401337792635,"#2A587F"],[0.36454849498327757,"#2A5880"],[0.36789297658862874,"#2A5980"],[0.37123745819397991,"#2A5981"],[0.37458193979933113,"#2A5982"],[0.3779264214046823,"#2B5A82"],[0.38127090301003352,"#2B5A83"],[0.38461538461538469,"#2B5B83"],[0.38795986622073592,"#2B5B84"],[0.39130434782608692,"#2C5C85"],[0.39464882943143809,"#2C5C85"],[0.39799331103678931,"#2C5D86"],[0.40133779264214048,"#2C5D86"],[0.40468227424749165,"#2C5D87"],[0.40802675585284287,"#2D5E87"],[0.41137123745819387,"#2D5E88"],[0.41471571906354504,"#2D5F89"],[0.41806020066889626,"#2D5F89"],[0.42140468227424743,"#2E608A"],[0.42474916387959866,"#2E608A"],[0.42809364548494983,"#2E618B"],[0.43143812709030105,"#2E618C"],[0.43478260869565222,"#2E618C"],[0.43812709030100339,"#2F628D"],[0.44147157190635461,"#2F628D"],[0.44481605351170578,"#2F638E"],[0.448160535117057,"#2F638F"],[0.451505016722408,"#30648F"],[0.45484949832775917,"#306490"],[0.4581939799331104,"#306590"],[0.46153846153846156,"#306591"],[0.46488294314381257,"#306592"],[0.46822742474916373,"#316692"],[0.47157190635451496,"#316693"],[0.47491638795986613,"#316793"],[0.47826086956521735,"#316794"],[0.48160535117056852,"#326895"],[0.48494983277591974,"#326895"],[0.48829431438127091,"#326996"],[0.49163879598662213,"#326996"],[0.4949832775919733,"#326997"],[0.49832775919732453,"#336A98"],[0.5016722408026757,"#336A98"],[0.50501672240802686,"#336B99"],[0.50836120401337792,"#336B99"],[0.51170568561872909,"#346C9A"],[0.51505016722408026,"#346C9B"],[0.51839464882943143,"#346D9B"],[0.52173913043478248,"#346D9C"],[0.52508361204013365,"#346E9D"],[0.52842809364548482,"#356E9D"],[0.5317725752508361,"#356E9E"],[0.53511705685618727,"#356F9E"],[0.53846153846153844,"#356F9F"],[0.5418060200668896,"#3670A0"],[0.54515050167224077,"#3670A0"],[0.54849498327759205,"#3671A1"],[0.55183946488294322,"#3671A1"],[0.55518394648829439,"#3772A2"],[0.55852842809364556,"#3772A3"],[0.56187290969899684,"#3773A3"],[0.56521739130434778,"#3773A4"],[0.56856187290969895,"#3773A4"],[0.57190635451505012,"#3874A5"],[0.5752508361204014,"#3874A6"],[0.57859531772575234,"#3875A6"],[0.58193979933110351,"#3875A7"],[0.58528428093645479,"#3976A8"],[0.58862876254180596,"#3976A8"],[0.59197324414715713,"#3977A9"],[0.5953177257525083,"#3977A9"],[0.59866220735785958,"#3978AA"],[0.60200668896321075,"#3A78AB"],[0.60535117056856191,"#3A79AB"],[0.60869565217391308,"#3A79AC"],[0.61204013377926425,"#3A79AC"],[0.61538461538461553,"#3B7AAD"],[0.6187290969899667,"#3B7AAE"],[0.62207357859531764,"#3B7BAE"],[0.62541806020066892,"#3B7BAF"],[0.62876254180602009,"#3C7CB0"],[0.63210702341137126,"#3C7CB0"],[0.63545150501672221,"#3C7DB1"],[0.63879598662207349,"#3C7DB1"],[0.64214046822742465,"#3C7EB2"],[0.64548494983277582,"#3D7EB3"],[0.64882943143812699,"#3D7FB3"],[0.65217391304347827,"#3D7FB4"],[0.65551839464882944,"#3D7FB5"],[0.65886287625418061,"#3E80B5"],[0.66220735785953178,"#3E80B6"],[0.66555183946488305,"#3E81B6"],[0.66889632107023422,"#3E81B7"],[0.67224080267558539,"#3F82B8"],[0.67558528428093656,"#3F82B8"],[0.67892976588628762,"#3F83B9"],[0.68227424749163879,"#3F83BA"],[0.68561872909698995,"#4084BA"],[0.68896321070234112,"#4084BB"],[0.69230769230769218,"#4085BB"],[0.69565217391304335,"#4085BC"],[0.69899665551839452,"#4086BD"],[0.70234113712374568,"#4186BD"],[0.70568561872909696,"#4186BE"],[0.70903010033444813,"#4187BF"],[0.7123745819397993,"#4187BF"],[0.71571906354515047,"#4288C0"],[0.71906354515050175,"#4288C1"],[0.72240802675585292,"#4289C1"],[0.72575250836120409,"#4289C2"],[0.72909698996655525,"#438AC2"],[0.73244147157190653,"#438AC3"],[0.7357859531772577,"#438BC4"],[0.73913043478260843,"#438BC4"],[0.74247491638796004,"#438CC5"],[0.74581939799331087,"#448CC6"],[0.74916387959866249,"#448DC6"],[0.75250836120401321,"#448DC7"],[0.75585284280936449,"#448EC8"],[0.75919732441471566,"#458EC8"],[0.76254180602006683,"#458FC9"],[0.76588628762541799,"#458FC9"],[0.76923076923076916,"#458FCA"],[0.77257525083612044,"#4690CB"],[0.77591973244147161,"#4690CB"],[0.77926421404682278,"#4691CC"],[0.78260869565217395,"#4691CD"],[0.78595317725752523,"#4792CD"],[0.78929765886287639,"#4792CE"],[0.79264214046822756,"#4793CF"],[0.7959866220735784,"#4793CF"],[0.79933110367893001,"#4894D0"],[0.80267558528428073,"#4894D0"],[0.80602006688963235,"#4895D1"],[0.80936454849498318,"#4895D2"],[0.81270903010033435,"#4896D2"],[0.81605351170568552,"#4996D3"],[0.81939799331103669,"#4997D4"],[0.82274247491638797,"#4997D4"],[0.82608695652173914,"#4998D5"],[0.8294314381270903,"#4A98D6"],[0.83277591973244147,"#4A99D6"],[0.83612040133779264,"#4A99D7"],[0.83946488294314392,"#4A9AD8"],[0.84280936454849509,"#4B9AD8"],[0.84615384615384626,"#4B9BD9"],[0.84949832775919742,"#4B9BDA"],[0.85284280936454826,"#4B9BDA"],[0.85618729096989987,"#4C9CDB"],[0.85953177257525071,"#4C9CDB"],[0.86287625418060221,"#4C9DDC"],[0.86622073578595304,"#4C9DDD"],[0.86956521739130421,"#4D9EDD"],[0.87290969899665538,"#4D9EDE"],[0.87625418060200666,"#4D9FDF"],[0.87959866220735783,"#4D9FDF"],[0.882943143812709,"#4DA0E0"],[0.88628762541806017,"#4EA0E1"],[0.88963210702341144,"#4EA1E1"],[0.89297658862876261,"#4EA1E2"],[0.89632107023411378,"#4EA2E3"],[0.89966555183946495,"#4FA2E3"],[0.90301003344481579,"#4FA3E4"],[0.9063545150501674,"#4FA3E5"],[0.90969899665551812,"#4FA4E5"],[0.91304347826086973,"#50A4E6"],[0.91638795986622057,"#50A5E7"],[0.91973244147157218,"#50A5E7"],[0.92307692307692291,"#50A6E8"],[0.92642140468227452,"#51A6E8"],[0.92976588628762535,"#51A7E9"],[0.93311036789297652,"#51A7EA"],[0.93645484949832769,"#51A8EA"],[0.93979933110367886,"#52A8EB"],[0.94314381270903014,"#52A9EC"],[0.94648829431438131,"#52A9EC"],[0.94983277591973247,"#52AAED"],[0.95317725752508364,"#53AAEE"],[0.95652173913043492,"#53ABEE"],[0.95986622073578565,"#53ABEF"],[0.96321070234113726,"#53ACF0"],[0.96655518394648809,"#54ACF0"],[0.96989966555183971,"#54ADF1"],[0.97324414715719043,"#54ADF2"],[0.97658862876254204,"#54AEF2"],[0.97993311036789288,"#55AEF3"],[0.98327759197324438,"#55AFF4"],[0.98662207357859522,"#55AFF4"],[0.98996655518394638,"#55B0F5"],[0.99331103678929766,"#56B0F6"],[0.99665551839464883,"#56B1F6"],[1,"#56B1F7"]],"colorbar":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"thickness":23.039999999999996,"title":"as.numeric(colData(se.subset)$RIN)","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":19.925280199252807},"tickmode":"array","ticktext":["5","6","7","8","9"],"tickvals":[0.068111111111111081,0.28959259259259257,0.51107407407407413,0.73255555555555563,0.95403703703703713],"tickfont":{"color":"rgba(0,0,0,1)","family":"","size":17.268576172685766},"ticklen":2,"len":0.5,"yanchor":"top","y":1}},"xaxis":"x","yaxis":"y","frame":null}],"layout":{"margin":{"t":23.305936073059364,"r":7.3059360730593621,"b":48.152760481527615,"l":76.048152760481543},"plot_bgcolor":"rgba(255,255,255,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.611872146118724},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.4859832355776485,0.57608061817640899],"tickmode":"array","ticktext":["-0.25","0.00","0.25","0.50"],"tickvals":[-0.25,0,0.25000000000000017,0.5],"categoryorder":"array","categoryarray":["-0.25","0.00","0.25","0.50"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":18.596928185969286},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"y","title":{"text":"PC1 (36.6%)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-0.52916375958825679,0.52780500537853048],"tickmode":"array","ticktext":["-0.50","-0.25","0.00","0.25","0.50"],"tickvals":[-0.5,-0.25,0,0.25,0.5],"categoryorder":"array","categoryarray":["-0.50","-0.25","0.00","0.25","0.50"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.6529680365296811,"tickwidth":0.66417600664176002,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":18.596928185969286},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":false,"gridcolor":null,"gridwidth":0,"zeroline":false,"anchor":"x","title":{"text":"PC2 (17.9%)","font":{"color":"rgba(0,0,0,1)","family":"","size":18.596928185969286}},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":"rgba(255,255,255,1)","line":{"color":"rgba(51,51,51,1)","width":0.66417600664176002,"linetype":"solid"},"yref":"paper","xref":"paper","layer":"below","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.8897637795275593,"font":{"color":"rgba(0,0,0,1)","family":"","size":17.268576172685766},"y":0.5,"yanchor":"top","title":{"text":"Number_Outliers_Binned<br />as.numeric(colData(se.subset)$RIN)","font":{"color":"rgba(0,0,0,1)","family":"","size":19.925280199252807}}},"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"source":"A","attrs":{"6920bd84582":{"x":{},"y":{},"colour":{},"label":{},"shape":{},"type":"scatter"}},"cur_data":"6920bd84582","visdat":{"6920bd84582":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.20000000000000001,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
   <td style="text-align:center;"> rnaseqc_Low.Quality.Reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 7 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> STAR_total_reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 8 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> STAR_uniquely_mapped </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 9 </td>
   <td style="text-align:center;"> PC1 </td>
   <td style="text-align:center;"> STAR_multimapped_multiple </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 10 </td>
   <td style="text-align:center;"> PC10 </td>
   <td style="text-align:center;"> heterozygosity_mean </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 11 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> rnaseqc_rRNA.Reads </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 12 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> rnaseqc_rRNA.Rate </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 13 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> rnaseqc_Median.3..bias </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 14 </td>
   <td style="text-align:center;"> PC2 </td>
   <td style="text-align:center;"> rnaseqc_Genes.Detected </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 15 </td>
   <td style="text-align:center;"> PC3 </td>
   <td style="text-align:center;"> rnaseqc_Median.Exon.CV </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 16 </td>
   <td style="text-align:center;"> PC3 </td>
   <td style="text-align:center;"> rnaseqc_Median.of.Avg.Transcript.Coverage </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 17 </td>
   <td style="text-align:center;"> PC3 </td>
   <td style="text-align:center;"> STAR_uniquely_mapped_percent </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 18 </td>
   <td style="text-align:center;"> PC3 </td>
   <td style="text-align:center;"> STAR_unmapped_tooshort_percent </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 19 </td>
   <td style="text-align:center;"> PC3 </td>
   <td style="text-align:center;"> DV200 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 20 </td>
   <td style="text-align:center;"> PC4 </td>
   <td style="text-align:center;"> RIN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 21 </td>
   <td style="text-align:center;"> PC5 </td>
   <td style="text-align:center;"> rnaseqc_Unique.Rate.of.Mapped </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 22 </td>
   <td style="text-align:center;"> PC5 </td>
   <td style="text-align:center;"> STAR_num_splices </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 23 </td>
   <td style="text-align:center;"> PC6 </td>
   <td style="text-align:center;"> STAR_multimapped_multiple_percent </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 24 </td>
   <td style="text-align:center;"> PC6 </td>
   <td style="text-align:center;"> STAR_avg_mapped_read_length </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 25 </td>
   <td style="text-align:center;"> PC7 </td>
   <td style="text-align:center;"> TIN_median </td>
  </tr>
  <tr>
   <td style="text-align:left;"> 26 </td>
   <td style="text-align:center;"> PC7 </td>
   <td style="text-align:center;"> rnaseqc_Duplicate.Rate.of.Mapped </td>
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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/heatmap_qc_outiers-1.png" style="display: block; margin: auto;" />

```r
out
```

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/heatmap_qc_outiers-2.png" style="display: block; margin: auto;" />

 

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
   <td style="text-align:left;"> SynCohort_P01_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P03_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P06_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P09_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P10_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P02_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P04_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P05_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P07_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P08_S1 </td>
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
   <td style="text-align:left;"> SynCohort_P01_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P03_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P06_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P09_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P10_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P02_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P07_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P08_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P04_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P05_S1 </td>
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
   <td style="text-align:left;"> SynCohort_P01_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P02_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P07_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P08_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P03_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P06_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P09_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P10_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P04_S1 </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P05_S1 </td>
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
   <td style="text-align:left;"> SynCohort_P01_S1 </td>
   <td style="text-align:center;"> 1 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P02_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P07_S1 </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P03_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P06_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P09_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P10_S1 </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P04_S1 </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P05_S1 </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P08_S1 </td>
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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/heatmap_sex_genes-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/sex_genes_pca-1.png" style="display: block; margin: auto;" />


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
   <td style="text-align:left;"> SynCohort_P01_S1 </td>
   <td style="text-align:center;"> SynCohort_P01 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 18359705 </td>
   <td style="text-align:center;"> 63.2 </td>
   <td style="text-align:center;"> 21497241 </td>
   <td style="text-align:center;"> 20165401 </td>
   <td style="text-align:center;"> 14332502 </td>
   <td style="text-align:center;"> 0.711 </td>
   <td style="text-align:center;"> 20033063 </td>
   <td style="text-align:center;"> 1464178 </td>
   <td style="text-align:center;"> 325372 </td>
   <td style="text-align:center;"> 0.0151 </td>
   <td style="text-align:center;"> 0.358 </td>
   <td style="text-align:center;"> 13319 </td>
   <td style="text-align:center;"> 0.842 </td>
   <td style="text-align:center;"> 0.225 </td>
   <td style="text-align:center;"> 24.41 </td>
   <td style="text-align:center;"> 21497241 </td>
   <td style="text-align:center;"> 19354063 </td>
   <td style="text-align:center;"> 90.03 </td>
   <td style="text-align:center;"> 1005386 </td>
   <td style="text-align:center;"> 4.68 </td>
   <td style="text-align:center;"> 1.71 </td>
   <td style="text-align:center;"> 142.2 </td>
   <td style="text-align:center;"> 2862762 </td>
   <td style="text-align:center;"> 9.0 </td>
   <td style="text-align:center;"> 66.7 </td>
   <td style="text-align:center;"> 0.605 </td>
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
   <td style="text-align:center;"> -4.8873349 </td>
   <td style="text-align:center;"> -0.6597682 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.471486 </td>
   <td style="text-align:center;"> 0.4806009 </td>
   <td style="text-align:center;"> 14.481787 </td>
   <td style="text-align:center;"> 6.885536 </td>
   <td style="text-align:center;"> 2.511932 </td>
   <td style="text-align:center;"> 2.4016849 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P02_S1 </td>
   <td style="text-align:center;"> SynCohort_P02 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 47132510 </td>
   <td style="text-align:center;"> 83.5 </td>
   <td style="text-align:center;"> 58941597 </td>
   <td style="text-align:center;"> 51412655 </td>
   <td style="text-align:center;"> 37794459 </td>
   <td style="text-align:center;"> 0.735 </td>
   <td style="text-align:center;"> 52801059 </td>
   <td style="text-align:center;"> 6140538 </td>
   <td style="text-align:center;"> 752332 </td>
   <td style="text-align:center;"> 0.0128 </td>
   <td style="text-align:center;"> 0.330 </td>
   <td style="text-align:center;"> 9123 </td>
   <td style="text-align:center;"> 0.950 </td>
   <td style="text-align:center;"> 0.428 </td>
   <td style="text-align:center;"> 19.75 </td>
   <td style="text-align:center;"> 58941597 </td>
   <td style="text-align:center;"> 49543142 </td>
   <td style="text-align:center;"> 84.05 </td>
   <td style="text-align:center;"> 2858833 </td>
   <td style="text-align:center;"> 4.85 </td>
   <td style="text-align:center;"> 2.89 </td>
   <td style="text-align:center;"> 102.3 </td>
   <td style="text-align:center;"> 1199569 </td>
   <td style="text-align:center;"> 9.2 </td>
   <td style="text-align:center;"> 89.2 </td>
   <td style="text-align:center;"> 0.739 </td>
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
   <td style="text-align:center;"> 3.0180342 </td>
   <td style="text-align:center;"> -3.1072062 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.612735 </td>
   <td style="text-align:center;"> 0.1700702 </td>
   <td style="text-align:center;"> 6.139679 </td>
   <td style="text-align:center;"> 12.516171 </td>
   <td style="text-align:center;"> 2.572288 </td>
   <td style="text-align:center;"> 1.1778128 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P03_S1 </td>
   <td style="text-align:center;"> SynCohort_P03 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 31813743 </td>
   <td style="text-align:center;"> 78.4 </td>
   <td style="text-align:center;"> 37270050 </td>
   <td style="text-align:center;"> 32944541 </td>
   <td style="text-align:center;"> 24843358 </td>
   <td style="text-align:center;"> 0.754 </td>
   <td style="text-align:center;"> 31296763 </td>
   <td style="text-align:center;"> 5973287 </td>
   <td style="text-align:center;"> 480734 </td>
   <td style="text-align:center;"> 0.0129 </td>
   <td style="text-align:center;"> 0.571 </td>
   <td style="text-align:center;"> 15432 </td>
   <td style="text-align:center;"> 0.636 </td>
   <td style="text-align:center;"> 0.223 </td>
   <td style="text-align:center;"> 7.37 </td>
   <td style="text-align:center;"> 37270050 </td>
   <td style="text-align:center;"> 31946321 </td>
   <td style="text-align:center;"> 85.72 </td>
   <td style="text-align:center;"> 720838 </td>
   <td style="text-align:center;"> 1.93 </td>
   <td style="text-align:center;"> 0.60 </td>
   <td style="text-align:center;"> 145.5 </td>
   <td style="text-align:center;"> 1749427 </td>
   <td style="text-align:center;"> 5.6 </td>
   <td style="text-align:center;"> 91.0 </td>
   <td style="text-align:center;"> 0.590 </td>
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
   <td style="text-align:center;"> -2.7767562 </td>
   <td style="text-align:center;"> 2.1486686 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.280046 </td>
   <td style="text-align:center;"> -0.2132391 </td>
   <td style="text-align:center;"> 13.913628 </td>
   <td style="text-align:center;"> 7.032218 </td>
   <td style="text-align:center;"> 2.320493 </td>
   <td style="text-align:center;"> 1.7395039 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P04_S1 </td>
   <td style="text-align:center;"> SynCohort_P04 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 50092444 </td>
   <td style="text-align:center;"> 82.2 </td>
   <td style="text-align:center;"> 58303064 </td>
   <td style="text-align:center;"> 53456545 </td>
   <td style="text-align:center;"> 43826329 </td>
   <td style="text-align:center;"> 0.820 </td>
   <td style="text-align:center;"> 47534295 </td>
   <td style="text-align:center;"> 10768769 </td>
   <td style="text-align:center;"> 298861 </td>
   <td style="text-align:center;"> 0.0051 </td>
   <td style="text-align:center;"> 0.574 </td>
   <td style="text-align:center;"> 15313 </td>
   <td style="text-align:center;"> 0.361 </td>
   <td style="text-align:center;"> 0.174 </td>
   <td style="text-align:center;"> 24.65 </td>
   <td style="text-align:center;"> 58303064 </td>
   <td style="text-align:center;"> 49993440 </td>
   <td style="text-align:center;"> 85.75 </td>
   <td style="text-align:center;"> 2272648 </td>
   <td style="text-align:center;"> 3.90 </td>
   <td style="text-align:center;"> 4.09 </td>
   <td style="text-align:center;"> 143.2 </td>
   <td style="text-align:center;"> 4930638 </td>
   <td style="text-align:center;"> 8.6 </td>
   <td style="text-align:center;"> 50.5 </td>
   <td style="text-align:center;"> 0.721 </td>
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
   <td style="text-align:center;"> 4.0530569 </td>
   <td style="text-align:center;"> 2.2036576 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.179169 </td>
   <td style="text-align:center;"> -0.2135360 </td>
   <td style="text-align:center;"> 7.508142 </td>
   <td style="text-align:center;"> 12.116255 </td>
   <td style="text-align:center;"> 2.138722 </td>
   <td style="text-align:center;"> 2.5681302 </td>
   <td style="text-align:center;"> Near boundary </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P05_S1 </td>
   <td style="text-align:center;"> SynCohort_P05 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 48793116 </td>
   <td style="text-align:center;"> 60.2 </td>
   <td style="text-align:center;"> 55510196 </td>
   <td style="text-align:center;"> 52059069 </td>
   <td style="text-align:center;"> 39010598 </td>
   <td style="text-align:center;"> 0.749 </td>
   <td style="text-align:center;"> 45216055 </td>
   <td style="text-align:center;"> 10294141 </td>
   <td style="text-align:center;"> 283922 </td>
   <td style="text-align:center;"> 0.0051 </td>
   <td style="text-align:center;"> 0.659 </td>
   <td style="text-align:center;"> 16851 </td>
   <td style="text-align:center;"> 0.751 </td>
   <td style="text-align:center;"> 0.119 </td>
   <td style="text-align:center;"> 7.48 </td>
   <td style="text-align:center;"> 55510196 </td>
   <td style="text-align:center;"> 41647977 </td>
   <td style="text-align:center;"> 75.03 </td>
   <td style="text-align:center;"> 2561539 </td>
   <td style="text-align:center;"> 4.61 </td>
   <td style="text-align:center;"> 1.00 </td>
   <td style="text-align:center;"> 98.2 </td>
   <td style="text-align:center;"> 2313096 </td>
   <td style="text-align:center;"> 7.5 </td>
   <td style="text-align:center;"> 67.3 </td>
   <td style="text-align:center;"> 0.632 </td>
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
   <td style="text-align:center;"> 2.3237907 </td>
   <td style="text-align:center;"> 3.1160063 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.885605 </td>
   <td style="text-align:center;"> -0.2999918 </td>
   <td style="text-align:center;"> 14.104713 </td>
   <td style="text-align:center;"> 6.232450 </td>
   <td style="text-align:center;"> 2.926052 </td>
   <td style="text-align:center;"> 2.7370001 </td>
   <td style="text-align:center;"> Within-sex outlier </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> TRUE </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P06_S1 </td>
   <td style="text-align:center;"> SynCohort_P06 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 42054109 </td>
   <td style="text-align:center;"> 75.6 </td>
   <td style="text-align:center;"> 45599151 </td>
   <td style="text-align:center;"> 44585310 </td>
   <td style="text-align:center;"> 39227893 </td>
   <td style="text-align:center;"> 0.880 </td>
   <td style="text-align:center;"> 38845330 </td>
   <td style="text-align:center;"> 6753821 </td>
   <td style="text-align:center;"> 382575 </td>
   <td style="text-align:center;"> 0.0084 </td>
   <td style="text-align:center;"> 0.597 </td>
   <td style="text-align:center;"> 13171 </td>
   <td style="text-align:center;"> 1.139 </td>
   <td style="text-align:center;"> 0.198 </td>
   <td style="text-align:center;"> 12.40 </td>
   <td style="text-align:center;"> 45599151 </td>
   <td style="text-align:center;"> 37442976 </td>
   <td style="text-align:center;"> 82.11 </td>
   <td style="text-align:center;"> 1556708 </td>
   <td style="text-align:center;"> 3.41 </td>
   <td style="text-align:center;"> 2.93 </td>
   <td style="text-align:center;"> 137.1 </td>
   <td style="text-align:center;"> 1683986 </td>
   <td style="text-align:center;"> 6.9 </td>
   <td style="text-align:center;"> 71.3 </td>
   <td style="text-align:center;"> 0.676 </td>
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
   <td style="text-align:center;"> 0.3457720 </td>
   <td style="text-align:center;"> 1.6524291 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.402021 </td>
   <td style="text-align:center;"> -0.1528213 </td>
   <td style="text-align:center;"> 7.130736 </td>
   <td style="text-align:center;"> 12.336318 </td>
   <td style="text-align:center;"> 2.361574 </td>
   <td style="text-align:center;"> 1.0647730 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 3 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P07_S1 </td>
   <td style="text-align:center;"> SynCohort_P07 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 52388122 </td>
   <td style="text-align:center;"> 40.2 </td>
   <td style="text-align:center;"> 58838664 </td>
   <td style="text-align:center;"> 55822638 </td>
   <td style="text-align:center;"> 39185883 </td>
   <td style="text-align:center;"> 0.702 </td>
   <td style="text-align:center;"> 53746914 </td>
   <td style="text-align:center;"> 5091750 </td>
   <td style="text-align:center;"> 1112443 </td>
   <td style="text-align:center;"> 0.0189 </td>
   <td style="text-align:center;"> 0.267 </td>
   <td style="text-align:center;"> 16519 </td>
   <td style="text-align:center;"> 0.875 </td>
   <td style="text-align:center;"> 0.240 </td>
   <td style="text-align:center;"> 24.24 </td>
   <td style="text-align:center;"> 58838664 </td>
   <td style="text-align:center;"> 51332417 </td>
   <td style="text-align:center;"> 87.24 </td>
   <td style="text-align:center;"> 2074668 </td>
   <td style="text-align:center;"> 3.53 </td>
   <td style="text-align:center;"> 3.07 </td>
   <td style="text-align:center;"> 117.2 </td>
   <td style="text-align:center;"> 2953020 </td>
   <td style="text-align:center;"> 8.1 </td>
   <td style="text-align:center;"> 87.5 </td>
   <td style="text-align:center;"> 0.699 </td>
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
   <td style="text-align:center;"> 2.6347340 </td>
   <td style="text-align:center;"> -2.6509361 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.481824 </td>
   <td style="text-align:center;"> -0.1436339 </td>
   <td style="text-align:center;"> 12.814842 </td>
   <td style="text-align:center;"> 6.624704 </td>
   <td style="text-align:center;"> 2.522271 </td>
   <td style="text-align:center;"> 0.2445507 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P08_S1 </td>
   <td style="text-align:center;"> SynCohort_P08 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 37710455 </td>
   <td style="text-align:center;"> 83.3 </td>
   <td style="text-align:center;"> 44753528 </td>
   <td style="text-align:center;"> 41336305 </td>
   <td style="text-align:center;"> 32815755 </td>
   <td style="text-align:center;"> 0.794 </td>
   <td style="text-align:center;"> 35804640 </td>
   <td style="text-align:center;"> 8948888 </td>
   <td style="text-align:center;"> 863276 </td>
   <td style="text-align:center;"> 0.0193 </td>
   <td style="text-align:center;"> 0.344 </td>
   <td style="text-align:center;"> 12428 </td>
   <td style="text-align:center;"> 0.653 </td>
   <td style="text-align:center;"> 0.164 </td>
   <td style="text-align:center;"> 21.87 </td>
   <td style="text-align:center;"> 44753528 </td>
   <td style="text-align:center;"> 40984763 </td>
   <td style="text-align:center;"> 91.58 </td>
   <td style="text-align:center;"> 2125588 </td>
   <td style="text-align:center;"> 4.75 </td>
   <td style="text-align:center;"> 3.29 </td>
   <td style="text-align:center;"> 98.1 </td>
   <td style="text-align:center;"> 1074750 </td>
   <td style="text-align:center;"> 4.7 </td>
   <td style="text-align:center;"> 51.5 </td>
   <td style="text-align:center;"> 0.680 </td>
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
   <td style="text-align:center;"> 0.7536962 </td>
   <td style="text-align:center;"> -1.3408014 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.624669 </td>
   <td style="text-align:center;"> 0.0032940 </td>
   <td style="text-align:center;"> 6.427962 </td>
   <td style="text-align:center;"> 12.558097 </td>
   <td style="text-align:center;"> 2.584223 </td>
   <td style="text-align:center;"> 1.9693726 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P09_S1 </td>
   <td style="text-align:center;"> SynCohort_P09 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 31193517 </td>
   <td style="text-align:center;"> 40.4 </td>
   <td style="text-align:center;"> 33337088 </td>
   <td style="text-align:center;"> 32018928 </td>
   <td style="text-align:center;"> 26530945 </td>
   <td style="text-align:center;"> 0.829 </td>
   <td style="text-align:center;"> 27851700 </td>
   <td style="text-align:center;"> 5485388 </td>
   <td style="text-align:center;"> 501965 </td>
   <td style="text-align:center;"> 0.0151 </td>
   <td style="text-align:center;"> 0.297 </td>
   <td style="text-align:center;"> 9579 </td>
   <td style="text-align:center;"> 0.847 </td>
   <td style="text-align:center;"> 0.222 </td>
   <td style="text-align:center;"> 10.58 </td>
   <td style="text-align:center;"> 33337088 </td>
   <td style="text-align:center;"> 27381231 </td>
   <td style="text-align:center;"> 82.13 </td>
   <td style="text-align:center;"> 1467476 </td>
   <td style="text-align:center;"> 4.40 </td>
   <td style="text-align:center;"> 3.72 </td>
   <td style="text-align:center;"> 143.1 </td>
   <td style="text-align:center;"> 2357940 </td>
   <td style="text-align:center;"> 7.6 </td>
   <td style="text-align:center;"> 91.5 </td>
   <td style="text-align:center;"> 0.707 </td>
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
   <td style="text-align:center;"> -2.1632964 </td>
   <td style="text-align:center;"> -0.9894507 </td>
   <td style="text-align:center;"> Female </td>
   <td style="text-align:center;"> -2.252013 </td>
   <td style="text-align:center;"> 0.2255427 </td>
   <td style="text-align:center;"> 13.911924 </td>
   <td style="text-align:center;"> 7.079181 </td>
   <td style="text-align:center;"> 2.292459 </td>
   <td style="text-align:center;"> 0.8772603 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 2 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> SynCohort_P10_S1 </td>
   <td style="text-align:center;"> SynCohort_P10 </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 28240502 </td>
   <td style="text-align:center;"> 50.8 </td>
   <td style="text-align:center;"> 33869930 </td>
   <td style="text-align:center;"> 29623711 </td>
   <td style="text-align:center;"> 20748229 </td>
   <td style="text-align:center;"> 0.700 </td>
   <td style="text-align:center;"> 32468250 </td>
   <td style="text-align:center;"> 1401680 </td>
   <td style="text-align:center;"> 505735 </td>
   <td style="text-align:center;"> 0.0149 </td>
   <td style="text-align:center;"> 0.592 </td>
   <td style="text-align:center;"> 12423 </td>
   <td style="text-align:center;"> 1.193 </td>
   <td style="text-align:center;"> 0.107 </td>
   <td style="text-align:center;"> 10.24 </td>
   <td style="text-align:center;"> 33869930 </td>
   <td style="text-align:center;"> 28184084 </td>
   <td style="text-align:center;"> 83.21 </td>
   <td style="text-align:center;"> 1124239 </td>
   <td style="text-align:center;"> 3.32 </td>
   <td style="text-align:center;"> 1.05 </td>
   <td style="text-align:center;"> 110.2 </td>
   <td style="text-align:center;"> 1118654 </td>
   <td style="text-align:center;"> 7.9 </td>
   <td style="text-align:center;"> 89.5 </td>
   <td style="text-align:center;"> 0.736 </td>
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
   <td style="text-align:center;"> -3.3016966 </td>
   <td style="text-align:center;"> -0.3725990 </td>
   <td style="text-align:center;"> Male </td>
   <td style="text-align:center;"> 2.552379 </td>
   <td style="text-align:center;"> 0.1437143 </td>
   <td style="text-align:center;"> 7.371112 </td>
   <td style="text-align:center;"> 12.587928 </td>
   <td style="text-align:center;"> 2.511932 </td>
   <td style="text-align:center;"> 1.2199114 </td>
   <td style="text-align:center;"> Pass </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> FALSE </td>
   <td style="text-align:center;"> 4 </td>
  </tr>
</tbody>
</table></div>

```r
qc_tbl %>%
  tibble::rownames_to_column("sample_id") %>%
  readr::write_tsv(qcFile)
```

### Heatmap of all QC flags (tissue/site subset)

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

<img src="/Users/ykefella/dev/htanBU-RNAseqQC/example_data/output/example_QCReport_files/figure-html/flagged_all_samples-1.png" style="display: block; margin: auto;" />

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
The filtered per tissue type SE object with appended RNA annotation will be saved to the file: ``Lung_SiteA_Gene_Expression.rds``.

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
##            used  (Mb) gc trigger  (Mb) limit (Mb) max used  (Mb)
## Ncells  8370215 447.1   16089522 859.3         NA 16089522 859.3
## Vcells 14879723 113.6   26728521 204.0      16384 20336076 155.2
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
##  [1] kableExtra_1.4.0            igraph_2.1.4               
##  [3] lubridate_1.9.2             forcats_1.0.0              
##  [5] stringr_1.5.0               dplyr_1.2.1                
##  [7] purrr_1.0.1                 readr_2.1.4                
##  [9] tidyr_1.3.2                 tibble_3.2.1               
## [11] tidyverse_2.0.0             plotly_4.12.0              
## [13] ggrepel_0.9.8               ComplexHeatmap_2.18.0      
## [15] ggpubr_0.6.1                readxl_1.4.2               
## [17] pheatmap_1.0.13             ggplot2_4.0.3              
## [19] circlize_0.4.18             RColorBrewer_1.1-3         
## [21] edgeR_4.0.16                limma_3.58.1               
## [23] rmarkdown_2.22              SummarizedExperiment_1.32.0
## [25] Biobase_2.62.0              GenomicRanges_1.54.1       
## [27] GenomeInfoDb_1.38.8         IRanges_2.36.0             
## [29] S4Vectors_0.40.2            BiocGenerics_0.48.1        
## [31] MatrixGenerics_1.14.0       matrixStats_1.5.0          
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
## [85] htmlwidgets_1.6.4       farver_2.1.1            htmltools_0.5.9        
## [88] lifecycle_1.0.5         httr_1.4.6              GlobalOptions_0.1.4    
## [91] statmod_1.5.0           bit64_4.0.5
```

