# =============================================================================
# Proteomics Mfuzz Processing
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(impute)
  library(Biobase)
  library(Mfuzz)
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
  annotation = file.path(root_dir, "01_Input", "00_OG1RF_DB.xlsx"),
  output = file.path(root_dir, "03_Results")
))

n_clusters <- 8
mfuzz_seed <- 123
transcript_mean_count_threshold <- 10

dir.create(paths$output, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# Helpers
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
# Proteomics Normalization
# -----------------------------------------------------------------------------

normalize_proteomics <- function(raw_mat) {
  # Retain proteins detected in at least 75% of samples.
  mat <- raw_mat[rowMeans(!is.na(raw_mat)) >= 0.75, , drop = FALSE]

  # Normalize each sample to total protein abundance and rescale to the mean sample total.
  sample_totals <- colSums(mat, na.rm = TRUE)
  mat <- sweep(mat, 2, sample_totals, "/") * mean(sample_totals, na.rm = TRUE)

  # Express each protein relative to its mean abundance across the time course, then log2-transform.
  protein_means <- rowMeans(mat, na.rm = TRUE)
  mat <- sweep(mat, 1, protein_means, "/")
  mat[mat <= 0] <- NA
  mat <- log2(mat)

  if (any(is.na(mat))) mat <- impute::impute.knn(mat)$data
  mat
}


# -----------------------------------------------------------------------------
# Mfuzz Clustering
# -----------------------------------------------------------------------------

run_mfuzz <- function(z_mat) {
  set.seed(mfuzz_seed)
  eset <- ExpressionSet(assayData = z_mat)
  eset <- standardise(eset)
  finite_row <- apply(exprs(eset), 1, function(x) all(is.finite(x)))
  eset <- eset[finite_row, ]

  m_val <- mestimate(eset)
  mfuzz_result <- mfuzz(eset, c = n_clusters, m = m_val)

  list(
    result = mfuzz_result,
    m = m_val
  )
}

make_cluster_membership <- function(mfuzz_result) {
  cluster <- as.integer(mfuzz_result$cluster)
  loci <- rownames(mfuzz_result$membership)

  data.frame(
    Locus_Tag = loci,
    Cluster = cluster,
    Membership = mfuzz_result$membership[cbind(loci, cluster)],
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Centroids
# -----------------------------------------------------------------------------

get_transcript_filter_loci <- function(input_path) {
  rna_raw <- read_counts_sheet(input_path, "Transcriptomics")
  dds <- DESeqDataSetFromMatrix(
    round(rna_raw),
    data.frame(dummy = factor(rep(1, ncol(rna_raw))), row.names = colnames(rna_raw)),
    design = ~ 1
  )
  dds <- estimateSizeFactors(dds)
  keep <- rowMeans(counts(dds, normalized = TRUE)) >= transcript_mean_count_threshold
  rownames(rna_raw)[keep]
}

calculate_centroids <- function(z_mat, cluster_df, loci_to_use) {
  centroids <- t(sapply(seq_len(n_clusters), function(k) {
    members <- cluster_df$Locus_Tag[cluster_df$Cluster == k]
    members <- intersect(members, loci_to_use)
    colMeans(z_mat[members, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(centroids) <- seq_len(n_clusters)
  centroids
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

prot_raw <- read_counts_sheet(paths$input, "Proteomics")
prot_hours <- extract_hours(colnames(prot_raw))

prot_norm <- normalize_proteomics(prot_raw)
prot_avg <- average_replicates_by_hour(prot_norm, prot_hours)
prot_z <- zscore_rows(prot_avg)

mfuzz_fit <- run_mfuzz(prot_z)
cluster_assignments <- make_cluster_membership(mfuzz_fit$result)

transcript_loci <- get_transcript_filter_loci(paths$input)
common_loci <- intersect(cluster_assignments$Locus_Tag, transcript_loci)
protein_centroids <- calculate_centroids(prot_z, cluster_assignments, common_loci)

annotation <- openxlsx::read.xlsx(paths$annotation) %>%
  mutate(Locus_Tag = as.character(Locus_Tag)) %>%
  select(Locus_Tag, Old_Locus_Tag, Gene_Name, Product)

cluster_members <- cluster_assignments %>%
  arrange(Cluster, desc(Membership)) %>%
  mutate(Membership = round(Membership, 4)) %>%
  left_join(annotation, by = "Locus_Tag")


# -----------------------------------------------------------------------------
# Write Results
# -----------------------------------------------------------------------------

out_file <- file.path(paths$output, "proteomics_mfuzz_results.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "raw")
writeData(wb, "raw", matrix_to_locus_df(prot_raw[rownames(prot_norm), , drop = FALSE]))
addWorksheet(wb, "normalized")
writeData(wb, "normalized", matrix_to_locus_df(prot_norm))
addWorksheet(wb, "normalized_avg")
writeData(wb, "normalized_avg", matrix_to_locus_df(prot_avg))
addWorksheet(wb, "z_score_mfuzz")
writeData(wb, "z_score_mfuzz", matrix_to_locus_df(prot_z))
addWorksheet(wb, "cluster_membership")
writeData(wb, "cluster_membership", cluster_members)
addWorksheet(wb, "protein_centroids")
writeData(wb, "protein_centroids", data.frame(Cluster = rownames(protein_centroids), protein_centroids, check.names = FALSE))
saveWorkbook(wb, out_file, overwrite = TRUE)

cat(sprintf(
  paste0(
    "Wrote %s\n",
    "Proteins retained after missingness filter: %d\n",
    "Mfuzz clusters: %d; m = %.3f; seed = %d\n"
  ),
  normalizePath(out_file), nrow(prot_z), n_clusters, mfuzz_fit$m, mfuzz_seed
))