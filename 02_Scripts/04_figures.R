# =============================================================================
# Figure Generation
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
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
  results = file.path(root_dir, "03_Results"),
  figures = file.path(root_dir, "05_Figures")
))

dir.create(paths$figures, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

read_centroid_sheet <- function(path, sheet) {
  if (!file.exists(path)) stop("Missing required results file: ", path)
  openxlsx::read.xlsx(path, sheet = sheet)
}

long_centroids <- function(df, modality) {
  df %>%
    mutate(Cluster = as.integer(Cluster)) %>%
    pivot_longer(-Cluster, names_to = "Time", values_to = "Z_score") %>%
    mutate(
      Time_h = as.numeric(sub("h$", "", Time)),
      Modality = modality
    )
}

plot_centroid_overlay <- function(plot_data) {
  ggplot(plot_data, aes(Time_h, Z_score, group = Modality, linetype = Modality, color = Modality)) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey78") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.7) +
    facet_wrap(~ Cluster, ncol = 4, labeller = label_both) +
    scale_x_continuous(breaks = sort(unique(plot_data$Time_h))) +
    scale_color_manual(values = c(Protein = "black", mRNA = "#D55E00")) +
    scale_linetype_manual(values = c(Protein = "solid", mRNA = "longdash")) +
    labs(x = "Time (h)", y = "Mean row z-score", color = NULL, linetype = NULL) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "top",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.8, "lines")
    )
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

proteomics_file <- file.path(paths$results, "proteomics_mfuzz_results.xlsx")
transcriptomics_file <- file.path(paths$results, "transcriptomics_mfuzz_results.xlsx")

protein_centroids <- read_centroid_sheet(proteomics_file, "protein_centroids")
mrna_centroids <- read_centroid_sheet(transcriptomics_file, "mrna_centroids")

plot_data <- bind_rows(
  long_centroids(protein_centroids, "Protein"),
  long_centroids(mrna_centroids, "mRNA")
) %>%
  mutate(
    Cluster = factor(Cluster, levels = sort(unique(Cluster))),
    Modality = factor(Modality, levels = c("Protein", "mRNA"))
  )

figure <- plot_centroid_overlay(plot_data)


# -----------------------------------------------------------------------------
# Write Results
# -----------------------------------------------------------------------------

write.csv(plot_data, file.path(paths$figures, "Figure_centroid_overlay_data.csv"), row.names = FALSE)

ggsave(
  file.path(paths$figures, "Figure_mfuzz_centroid_overlay.pdf"),
  figure,
  width = 7.2,
  height = 4.8,
  units = "in"
)

ggsave(
  file.path(paths$figures, "Figure_mfuzz_centroid_overlay.png"),
  figure,
  width = 7.2,
  height = 4.8,
  units = "in",
  dpi = 300
)

cat(sprintf("Wrote Figures to %s\n", normalizePath(paths$figures)))
