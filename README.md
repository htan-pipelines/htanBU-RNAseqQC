# htanBU-RNAseqQC

A downstream, cohort-level QC pipeline for bulk RNA-seq data processed with the
[htan-pipelines/bulk-rna-seq-pipeline](https://github.com/htan-pipelines/bulk-rna-seq-pipeline).
It picks up where that pipeline leaves off: given the aggregated
`SummarizedExperiment` (plus a few sidecar files), it produces an HTML QC report
covering sample identity/contamination checks, per-metric QC scatterplots,
expression-based outlier detection, QC-metric clustering, a sex-check from
X/Y-linked gene expression, and a combined per-sample QC flag summary — plus a
standalone interactive sample-relatedness network.

## Where this fits in the pipeline

```
per-sample:  RNA_seq_pipeline.wdl   (FastQC, STAR, RSEM, RNA-SeQC2, RSeQC TIN,
                                      samtools, somalier extract, arcasHLA genotype)
                    |
cohort-wide: Aggregation/aggregate.wdl (somalier_relate/final, arcasHLA_merge,
                                          combine_se  ->  aggregated SummarizedExperiment,
                                          somalier.pairs.tsv, genotypes.tsv)
                    |
             *** htanBU-RNAseqQC (this repo) ***
                    |
             qc_report.Rmd (HTML)  +  somalier_network_graph.html (interactive network)
```

Run both `RNA_seq_pipeline.wdl` (once per sample) and `Aggregation/aggregate.wdl`
(once for the whole cohort) first. This repo does not touch raw reads, BAMs, or
per-sample outputs directly — everything it needs comes out of the aggregation step.

## Required inputs

| Param | What it is | Produced by | How to load it |
|---|---|---|---|
| `se` | Aggregated `SummarizedExperiment`. Assays: `expected_count`, `TPM`, `FPKM`. `rowData` has 8 columns including `gene_id`, `hgnc_symbol`, `entrezgene_id`, `transcript_count`. `colData` has one row per sample with QC columns prefixed `fastqc_`, `STAR_`, `rnaseqc_`, `samtools_`, `TIN_`, plus `Somalier`/`Genotypes` summary columns. | `Aggregation/aggregate.wdl` → `combine_se.wdl` | `se <- readRDS("Aggregated_Gene_Expression.rds")` |
| `rnaAnnot` | Your own per-sample metadata table: sample ID, patient/participant ID, tissue type, collection site, cohort, sequencing batch, RIN, DV200, etc. This is **not** produced by the pipeline — it's whatever sample-tracking sheet you already keep (e.g. a Terra data table export). | You / your LIMS / Terra | `rnaAnnot <- read.table("metadata.tsv", sep = "\t", header = TRUE)` |
| `somalier_pairs` | Pairwise sample relatedness, used for sample-swap/contamination and identity checks. | `Aggregation/aggregate.wdl` → `somalier_final.wdl` | `somalier_pairs <- read.delim("somalier.pairs.tsv")` |
| `genotypes` | Per-sample HLA genotype calls, used as an independent identity check. | `Aggregation/aggregate.wdl` → `arcasHLA_merge.wdl` | `genotypes <- read.delim("genotypes.tsv")` |

That's it — four objects. Everything else (QC cutoffs, which metadata columns
mean what, whether to subset by tissue, whether to highlight a priority
sample subset) is a `param` with a sensible default, documented in the YAML
header of `qc_report.Rmd`.

### Column mapping

Your `rnaAnnot` file almost certainly doesn't use the exact column names this
report expects internally. Instead of renaming your file, set `params$column_map`
(a named list) to point at whatever your columns are actually called, e.g.:

```r
column_map = list(
  sample_id       = "your_sample_id_column",
  patient_id      = "your_patient_column",
  tissue_type     = "your_tissue_column",
  collection_site = "your_site_column",
  cohort          = "your_cohort_column",
  batch_id        = "your_batch_column",
  rin             = "rin",
  dv200           = "dv200"
)
```

Only `sample_id` is strictly required, and it must match `colnames(se)` after
mapping. Leave any other entry pointing at a column you don't have — the
corresponding plots will just show blank/NA groupings rather than failing.

## Repo layout

```
htanBU-RNAseqQC/
├── README.md
├── qc_report.Rmd                  # the QC report itself (parameterized R Markdown)
├── render_qc_report.R             # example driver: reads inputs, renders per tissue
├── generate_somalier_network.R    # standalone visNetwork -> somalier_network_graph.html
└── example_data/
    └── build_example_and_render.R # builds a synthetic SE + metadata/somalier/genotypes
                                    # tables and renders qc_report.Rmd against them end to
                                    # end -- a smoke test with no real data required
```

## Usage

Don't have real pipeline outputs handy yet? `Rscript example_data/build_example_and_render.R`
builds a small synthetic dataset (fake SE, metadata, somalier pairs, genotypes) and
renders the full report end to end against it, writing everything to
`example_data/output/` -- useful for confirming your R environment is set up
correctly before pointing the pipeline at real data.

1. Edit the `config` block at the top of `render_qc_report.R` to point at your
   four input files and (if needed) your `column_map`.
2. `Rscript render_qc_report.R`

This renders one HTML report per `tissue_type` found in your metadata, writes
the per-sample QC flag table as a `.tsv`, optionally writes the tissue-subset
`SummarizedExperiment`, and generates
`somalier_network_graph.html` — an interactive relatedness network you can open directly in
a browser (search/filter by sample or patient, hover for relatedness values,
cross-patient matches are colored red).

You can also render a single report by hand from R/RStudio:

```r
rmarkdown::render(
  "qc_report.Rmd",
  params = list(
    se = se, rnaAnnot = rnaAnnot, somalier_pairs = somalier_pairs, genotypes = genotypes,
    column_map = list(sample_id = "your_id_col", patient_id = "your_patient_col",
                       tissue_type = "your_tissue_col", collection_site = "your_site_col",
                       cohort = "your_cohort_col", batch_id = "your_batch_col"),
    tissue = "Your Tissue Type",
    tissueType = "My Cohort QC"
  )
)
```

### R package requirements

CRAN: `tidyverse`, `RColorBrewer`, `circlize`, `ggplot2`, `pheatmap`, `ggpubr`,
`ggrepel`, `plotly`, `igraph`, `kableExtra`, `readxl`, `Matrix`, `rmarkdown`,
`knitr`, `plyr` (used directly for `plyr::ldply()` in the QC-metrics-vs-PC-components
table; safe alongside `dplyr` since it's only ever called via `::`, never attached).
Bioconductor: `SummarizedExperiment`, `edgeR`, `ComplexHeatmap`. `edgeR` (via `limma`)
additionally requires the CRAN package `statmod`, which isn't always pulled in
automatically depending on what's already installed -- install it explicitly if
`library(edgeR)` fails with "there is no package called 'statmod'". For
`generate_somalier_network.R`: `dplyr` (part of tidyverse), `visNetwork`,
`htmlwidgets`.

Rendering to HTML also requires [pandoc](https://pandoc.org/installing.html)
(bundled with RStudio; install separately if running via plain `Rscript`).

```r
install.packages(c("tidyverse", "RColorBrewer", "circlize", "pheatmap", "ggpubr",
                    "ggrepel", "plotly", "igraph", "kableExtra", "readxl", "Matrix",
                    "rmarkdown", "knitr", "plyr", "statmod", "visNetwork", "htmlwidgets"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "edgeR", "ComplexHeatmap"))
```

## QC cutoffs used in this report

Set via `params$qc_cutoffs`; defaults shown below (all overridable, no need to
edit `qc_report.Rmd`):

| Cutoff | Default | Flags a sample when... |
|---|---|---|
| `tin_median` / `rin` | 50 / 5 | TIN < 50 **and** RIN < 5 (both low) |
| `threeprime_bias` | 0.5 | RNA-SeQC median 3' bias > 0.5 |
| `genes_detected` / `exon_cv` | 5000 / 1.0 | Genes detected < 5000 **or** median exon CV > 1.0 |
| `rrna_rate` | 0.01 | rRNA rate > 1% |
| `dv200` | 50 | DV200 < 50 |
| `heterozygosity_mean` | 0.5 | Somalier mean heterozygosity rate > 0.5 |
| `somalier_relatedness` | 0.6 | Cross-patient somalier relatedness ≥ 0.6 (possible swap) |
| `hla_match` | 0.6 | Cross-patient HLA genotype match ≥ 0.6 (possible swap) |

