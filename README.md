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
mean what, whether to subset by tissue/site, whether to highlight a priority
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
├── render_qc_report.R             # example driver: reads inputs, renders per tissue/site
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

This renders one HTML report per `tissue_type` × `collection_site` combination
found in your metadata (set `site_col = NULL` in the config to render one
report per tissue type only), writes the per-sample QC flag table as a `.tsv`,
optionally writes the tissue/site-subset `SummarizedExperiment`, and generates
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
    tissue = "Your Tissue Type", site = NULL,
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

## Changes made while generalizing this from the original PCGA02 script

The original script (`PCGA02_BulkRNA_QC_v2.Rmd` / `render_PCGA02_BulkRNA_QC_v2.R`)
was written for one specific cohort and had accumulated some redundancy and a
few real bugs, which are called out here explicitly rather than silently
fixed, per the request to make redundancy visible:

**Removed hardcoded, environment-specific paths** — `setwd()`, `knitr::opts_knit$set(root.dir = ...)`,
and all `/restricted/projectnb/...` absolute paths. The report and driver
script now take all paths as params/config.

**Deduplicated library loading** — `RColorBrewer`, `ggpubr`, and `readxl` were
each `library()`-ed twice.

**Trimmed unused dependencies** — of the 33 packages the original script
loaded, 16 (`DESeq2`, `sva`, `patchwork`, `gridExtra`, `lubridate`, `stringr`,
`stringi`, `viridis`, `fgsea`, `ggprism`, `rstatix`, `glue`, `Hmisc`,
`mixtools`, `plotmm`, `fields`, `ggraph`, `jsonlite`) were never actually
called anywhere in the report — almost certainly carried over from a shared
lab template. Dropping them cuts install time and dependency risk
significantly for a pipeline meant for public use (several, like `DESeq2`,
`sva`, and `fgsea`, are heavy Bioconductor installs unrelated to what this
report does).

**Removed ~130 lines of fully commented-out dead code** — an HLA
genotype-vs-JSON-files cross-check block that was entirely commented out and
superseded by the working `genotypes.tsv`-based version later in the same
report.

**Fixed: `genotypes` param silently ignored.** The HLA-comparison chunk
re-read a hardcoded `"TerraOutput/genotypes.tsv"` path instead of using
`params$genotypes` that had already been passed in — meaning the report used
whatever happened to be on disk at that relative path, not the file you told
it to use. Now uses `params$genotypes` throughout.

**Fixed: mislabeled sex-marker gene annotation.** The sex-gene heatmap
hardcoded `Gender <- factor(c(rep("Male", 6), "Female"))` — i.e. it assumed
exactly 6 Y-linked marker genes always survive filtering, in that exact
order, followed by XIST. If any Y-linked gene was missing from the SE
(entirely plausible on a different reference/annotation), the row labels
would silently shift and mislabel genes. Row labels are now derived from
which marker genes are actually present.

**Fixed: `output_se` was documented but never written.** The "save full SE"
chunk had its `saveRDS()` call commented out, so the report text promised a
file that was never created. Restored to actually respect `params$output_se`
(matching the already-working per-tissue save).

**Replaced a magic number.** The HLA genotype match score divided shared
fields by a hardcoded `14` (assumed HLA panel width); now uses `ncol(x)` so
it adapts to whatever HLA-typing panel width your genotypes.tsv has.

**Replaced fixed-length hardcoded color palettes** (hex vectors sized for
exactly 4 PCGA tissue types / 5 sites / 5 cohorts, redefined inconsistently
across sections) with one dynamic palette generator (`mk_disc_cols()`),
defined once and reused, that always produces exactly as many colors as
there are groups in your data.

**Made the "priority sample subset" heatmap fully optional.** It was
hardcoded to read a PCGA/UCL-specific Excel file
(`QC/PCABulkDataIDs_260114.xlsx`) with fixed column names
(`"RNA biospecimen ID"`, `"Parent biospecimen ID"`, `"Priority"`). This
section is now skipped entirely unless you set `params$priority_list_file`,
with the relevant column names configurable via `priority_id_cols`,
`priority_column`, and `priority_value`.

**Replaced the interactive somalier network's embedding.** The Rmd
previously drew two *static* igraph network plots (by patient ID and by TIN)
duplicating what an interactive network shows better. Per the request to
output this as a separate deliverable, the network is now built once by
`generate_somalier_network.R` and saved as a standalone `somalier_network_graph.html` —
not embedded in the knitted report — so it can be opened, searched, and
filtered on its own.

**Replaced the hardcoded, 4-tissue-type render loop.** The original driver
script had a manual `if/else` chain matching directory names to exactly the
four PCGA tissue types (`Bronchial Biopsy`, `Bronchial Brush`, `Nasal Brush`,
`Resection`). `render_qc_report.R` now derives the tissue × site
combinations to render directly from your metadata, so it works for any
cohort without editing the driver script.

**Removed the one-off sample-ID-correction hack** (a specific
`PCGA02_10044_1001842` → `PCGA02_10044_2101842` fix baked into the render
script) from the main path. Real cohorts do occasionally need this kind of
fix, so the pattern is preserved as a clearly-labeled, commented-out example
in `render_qc_report.R` rather than silently running by default.

**Generalized all metadata column references** (`entity.pcga2_biospecimen_id`,
`patient_ID`, `sequencing_batch_id`, etc.) to canonical names (`sample_id`,
`patient_id`, `batch_id`, ...) applied once via `params$column_map`, so any
lab's metadata sheet works without editing `qc_report.Rmd`.

## Publishing to `htan-pipelines/htanBU-RNAseqQC`

Suggested next steps once you're happy with the contents of this folder:

```bash
cd htanBU-RNAseqQC
git init
git add .
git commit -m "Initial generalized htanBU-RNAseqQC pipeline"
git remote add origin git@github.com:htan-pipelines/htanBU-RNAseqQC.git
git push -u origin main
```

Consider also adding: a small `example_data/` with a synthetic/subsetted SE +
metadata so new users can smoke-test the report without real data, a GitHub
Actions workflow that renders the example on every PR, and a `DESCRIPTION`
or `renv.lock` to pin package versions.
