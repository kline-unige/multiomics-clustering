# =============================================================================
# Transcriptomics Mfuzz Processing
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(DESeq2)
})


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else normalizePath(getwd())
root_dir <- normalizePath(file.path(script_dir, ".."))

parse_cli <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  for (arg in args) {
    if (!grepl("^--", arg)) next
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    if (length(kv) == 2) out[[kv[[1]]]] <- kv[[2]]
  }
  out
}

paths <- parse_cli(list(
  input = file.path(root_dir, "01_Input", "01_raw_counts_data.xlsx"),
  output = file.path(root_dir, "03_Results")
))

n_clusters <- 8
mean_count_threshold <- 10

dir.create(paths$output, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# General helpers
# -----------------------------------------------------------------------------

extract_hours <- function(sample_names) {
  match_pos <- regexpr("([0-9]+)h", sample_names, perl = TRUE)
  hours <- rep(NA_real_, length(sample_names))
  has_hour <- match_pos > 0
  if (any(has_hour)) {
    hours[has_hour] <- as.numeric(sub("h", "", regmatches(sample_names, match_pos)))
  }
  hours
}

read_counts_sheet <- function(path, sheet) {
  df <- openxlsx::read.xlsx(path, sheet = sheet, na.strings = c("", "NA"))

  id_col <- which(tolower(names(df)) %in% c("gene", "locus_tag"))
  if (!length(id_col)) id_col <- 1
  names(df)[id_col] <- "Locus_Tag"

  df <- df %>%
    mutate(Locus_Tag = trimws(sub("^\ufeff", "", Locus_Tag))) %>%
    filter(!is.na(Locus_Tag) & Locus_Tag != "") %>%
    filter(if_any(-Locus_Tag, ~ !is.na(.x))) %>%
    mutate(across(-Locus_Tag, ~ suppressWarnings(as.numeric(.x))))

  out <- as.data.frame(df)
  rownames(out) <- out$Locus_Tag
  as.matrix(out[, -1, drop = FALSE])
}

average_replicates_by_hour <- function(mat, hours) {
  groups <- split(seq_along(hours), hours)
  out <- sapply(groups, function(idx) {
    rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
  })

  ord <- order(as.numeric(names(groups)))
  out <- out[, ord, drop = FALSE]
  colnames(out) <- paste0(as.numeric(names(groups))[ord], "h")
  out
}

zscore_rows <- function(mat) {
  t(scale(t(as.matrix(mat)), center = TRUE, scale = TRUE))
}

matrix_to_locus_df <- function(mat) {
  data.frame(Locus_Tag = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE)
}


# -----------------------------------------------------------------------------
# Transcriptomics Normalization
# -----------------------------------------------------------------------------

normalize_transcriptomics <- function(raw_counts) {
  # Estimate DESeq2 size factors on all genes, then use normalized counts to apply low-expression filter.
  dds_all <- DESeqDataSetFromMatrix(
    round(raw_counts),
    data.frame(dummy = factor(rep(1, ncol(raw_counts))), row.names = colnames(raw_counts)),
    design = ~ 1
  )
  dds_all <- estimateSizeFactors(dds_all)
  keep <- rowMeans(counts(dds_all, normalized = TRUE)) >= mean_count_threshold

  # Variance-stabilizing transformation.
  dds_filtered <- DESeqDataSetFromMatrix(
    round(raw_counts[keep, , drop = FALSE]),
    data.frame(row.names = colnames(raw_counts)),
    design = ~ 1
  )

  list(
    normalized = assay(varianceStabilizingTransformation(dds_filtered, blind = TRUE)),
    retained_loci = rownames(raw_counts)[keep]
  )
}

calculate_centroids <- function(z_mat, cluster_df) {
  common_loci <- intersect(cluster_df$Locus_Tag, rownames(z_mat))

  centroids <- t(sapply(seq_len(n_clusters), function(k) {
    members <- cluster_df$Locus_Tag[cluster_df$Cluster == k]
    members <- intersect(members, common_loci)
    colMeans(z_mat[members, , drop = FALSE], na.rm = TRUE)
  }))

  rownames(centroids) <- seq_len(n_clusters)
  centroids
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

rna_raw <- read_counts_sheet(paths$input, "Transcriptomics")
rna_hours <- extract_hours(colnames(rna_raw))

rna_processed <- normalize_transcriptomics(rna_raw)
rna_norm <- rna_processed$normalized
rna_avg <- average_replicates_by_hour(rna_norm, rna_hours)
rna_z <- zscore_rows(rna_avg)

proteomics_results <- file.path(paths$output, "proteomics_mfuzz_results.xlsx")
if (!file.exists(proteomics_results)) {
  stop("Missing proteomics results.")
}

cluster_members <- openxlsx::read.xlsx(proteomics_results, sheet = "cluster_membership") %>%
  mutate(Cluster = as.integer(Cluster))

mrna_centroids <- calculate_centroids(rna_z, cluster_members)


# -----------------------------------------------------------------------------
# Write Results
# -----------------------------------------------------------------------------

out_file <- file.path(paths$output, "transcriptomics_mfuzz_results.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "raw")
writeData(wb, "raw", matrix_to_locus_df(rna_raw))
addWorksheet(wb, "normalized")
writeData(wb, "normalized", matrix_to_locus_df(rna_norm))
addWorksheet(wb, "normalized_avg")
writeData(wb, "normalized_avg", matrix_to_locus_df(rna_avg))
addWorksheet(wb, "z_score_mfuzz")
writeData(wb, "z_score_mfuzz", matrix_to_locus_df(rna_z))
addWorksheet(wb, "mrna_centroids")
writeData(wb, "mrna_centroids", data.frame(Cluster = rownames(mrna_centroids), mrna_centroids, check.names = FALSE))
saveWorkbook(wb, out_file, overwrite = TRUE)

cat(sprintf(
  paste0(
    "Wrote %s\n",
    "Transcripts retained after mean normalized count filter >= %s: %d\n"
  ),
  normalizePath(out_file), mean_count_threshold, nrow(rna_z)
))