#!/usr/bin/env Rscript
# ==============================================================================
# generate_somalier_network.R
#
# Builds an interactive sample-relatedness network from somalier pairwise output
# (Aggregation/aggregate.wdl -> somalier_final.wdl, "somalier.pairs.tsv") and
# writes it to a STANDALONE HTML file via visNetwork/htmlwidgets.
#
# This is deliberately kept OUT of qc_report.Rmd: htmlwidgets networks don't
# degrade gracefully inside a long, code-folded, multi-tab Rmd report the way a
# static plot does, and users often want to share/open just the network view on
# its own (e.g. to click through and find sample swaps).
#
# Usage as a script:
#   Rscript generate_somalier_network.R somalier.pairs.tsv somalier_network_graph.html [cutoff]
#
# Usage as a library (e.g. from render_qc_report.R):
#   source("generate_somalier_network.R")
#   generate_somalier_network(somalier_pairs, "somalier_network_graph.html", cutoff = 0.6)
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(visNetwork)
  library(htmlwidgets)
})

#' Default sample -> patient ID heuristic used by the original PCGA report:
#' takes the first two underscore-delimited tokens of the sample ID
#' (e.g. "PCGA02_10044_1001842" -> "PCGA02_10044"). This is a naming-convention
#' assumption, NOT something the bulk-rna-seq-pipeline itself defines -- override
#' `sample_to_patient` with your own function (or a named lookup vector wrapped
#' in a function) if your sample IDs don't follow a "prefix_patient_specimen"
#' pattern.
default_sample_to_patient <- function(sample_ids) {
  sapply(strsplit(sample_ids, "_"), function(x) paste(x[seq_len(min(2, length(x)))], collapse = "_"))
}

#' Generate an interactive somalier relatedness network and save it as a
#' standalone HTML file.
#'
#' @param somalier_pairs data.frame read from somalier.pairs.tsv (must have
#'   columns X.sample_a / sample_a, sample_b, relatedness -- read.delim() will
#'   turn a leading "#sample_a" header into "X.sample_a" automatically).
#' @param output_html Path to write the standalone HTML widget to.
#' @param cutoff Minimum relatedness to draw an edge for (default 0.6). Should
#'   match params$qc_cutoffs$somalier_relatedness in qc_report.Rmd if you want
#'   the two to agree.
#' @param sample_to_patient Function mapping a character vector of sample IDs
#'   to patient/participant IDs, used only for node coloring/grouping and to
#'   flag cross-patient edges. Defaults to `default_sample_to_patient()`.
generate_somalier_network <- function(somalier_pairs,
                                       output_html = "somalier_network_graph.html",
                                       cutoff = 0.6,
                                       sample_to_patient = default_sample_to_patient) {

  sample_a_col <- if ("X.sample_a" %in% names(somalier_pairs)) "X.sample_a" else "sample_a"
  stopifnot(sample_a_col %in% names(somalier_pairs), "sample_b" %in% names(somalier_pairs),
            "relatedness" %in% names(somalier_pairs))

  samples <- unique(c(somalier_pairs[[sample_a_col]], somalier_pairs$sample_b))
  patient_ids <- sample_to_patient(samples)
  names(patient_ids) <- samples
  patient_levels <- sort(unique(patient_ids))
  pal <- grDevices::hcl.colors(length(patient_levels), palette = "Set3")
  patient_color_map <- setNames(pal, patient_levels)

  nodes <- data.frame(
    id    = samples,
    label = "",  # keep blank to prevent overlap; use hover tooltip instead
    group = patient_ids[samples],
    color = patient_color_map[patient_ids[samples]],
    title = sprintf(
      "<b>Sample:</b> %s<br/><b>Patient:</b> %s",
      samples, patient_ids[samples]
    ),
    stringsAsFactors = FALSE
  )

  edges_df <- somalier_pairs %>%
    transmute(from = .data[[sample_a_col]], to = sample_b, rel = relatedness) %>%
    filter(!is.na(rel), rel >= cutoff) %>%
    filter(from %in% nodes$id, to %in% nodes$id) %>%
    distinct(from, to, .keep_all = TRUE)

  rel_min <- cutoff
  rel_max <- if (nrow(edges_df) > 0) max(edges_df$rel, na.rm = TRUE) else rel_min
  if (!is.finite(rel_max) || rel_max == rel_min) rel_max <- rel_min + 1e-6
  scale_width <- function(x) 1 + 7 * (x - rel_min) / (rel_max - rel_min)

  cross_patient <- patient_ids[edges_df$from] != patient_ids[edges_df$to]

  edges <- data.frame(
    from  = edges_df$from,
    to    = edges_df$to,
    width = scale_width(edges_df$rel),
    color = ifelse(cross_patient, "red", "rgba(0,0,0,0.25)"),
    title = sprintf(
      "<b>Relatedness:</b> %.3f<br/><b>From patient:</b> %s<br/><b>To patient:</b> %s",
      edges_df$rel,
      patient_ids[edges_df$from],
      patient_ids[edges_df$to]
    ),
    # NOTE: a bare `smooth = TRUE` here breaks data.frame() whenever no pairs clear
    # `cutoff` (edges_df has 0 rows): every other column is then a 0-length vector,
    # and a length-1 scalar can't be recycled down to 0 rows ("arguments imply
    # differing number of rows: 0, 1"). A cohort with no cross-sample relatedness
    # above cutoff is a real, unremarkable case, so this must not error.
    smooth = rep(TRUE, nrow(edges_df)),
    stringsAsFactors = FALSE
  )

  network <- visNetwork(nodes, edges, height = "800px", width = "100%") %>%
    visOptions(
      highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
      nodesIdSelection = list(enabled = TRUE, useLabels = FALSE, main = "Find sample"),
      selectedBy = list(variable = "group", main = "Filter by patient")
    ) %>%
    visInteraction(
      hover = TRUE,
      navigationButtons = TRUE,
      zoomView = TRUE,
      dragView = TRUE
    ) %>%
    visPhysics(
      solver = "forceAtlas2Based",
      stabilization = list(enabled = TRUE, iterations = 2000),
      forceAtlas2Based = list(gravitationalConstant = -30, springLength = 120, springConstant = 0.02)
    ) %>%
    visLayout(randomSeed = 123)

  htmlwidgets::saveWidget(network, file = output_html, selfcontained = TRUE)
  message("Wrote somalier relatedness network to: ", normalizePath(output_html, mustWork = FALSE))
  invisible(network)
}

# ------------------------------------------------------------------------------
# CLI entrypoint: only runs when this file is executed directly with Rscript,
# not when source()'d from another script (e.g. render_qc_report.R).
# ------------------------------------------------------------------------------
if (sys.nframe() == 0 && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript generate_somalier_network.R <somalier.pairs.tsv> <output.html> [cutoff]")
  }
  pairs_path <- args[1]
  output_html <- args[2]
  cutoff <- if (length(args) >= 3) as.numeric(args[3]) else 0.6

  somalier_pairs <- read.delim(pairs_path, header = TRUE, stringsAsFactors = FALSE)
  generate_somalier_network(somalier_pairs, output_html = output_html, cutoff = cutoff)
}
