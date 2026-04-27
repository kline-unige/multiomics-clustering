# =============================================================================
# Pathway Analysis
# =============================================================================

suppressPackageStartupMessages({
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
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
  annotation = file.path(root_dir, "01_Input", "00_OG1RF_DB.xlsx"),
  output = file.path(root_dir, "03_Results")
))

n_clusters <- 8
alpha_core <- 0.7
kegg_organism <- "efi"

dir.create(paths$output, showWarnings = FALSE, recursive = TRUE)


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

ratio_numerator <- function(x) as.numeric(sub("/.*$", "", x))
ratio_denominator <- function(x) as.numeric(sub("^.*/", "", x))

extract_go_terms <- function(annotation, id_col, name_col) {
  annotation %>%
    select(Locus_Tag, !!sym(id_col), !!sym(name_col)) %>%
    filter(!is.na(!!sym(id_col))) %>%
    mutate(
      GO_ID = str_split(!!sym(id_col), " ?// ?"),
      GO_Name = str_split(!!sym(name_col), " ?// ?")
    ) %>%
    tidyr::unnest(c(GO_ID, GO_Name)) %>%
    filter(GO_ID != "") %>%
    distinct() %>%
    transmute(Gene = Locus_Tag, GO_ID, GO_Name)
}

run_go_enrichment <- function(genes, term2gene, term2name, universe) {
  if (!length(genes)) return(data.frame())

  result <- tryCatch(
    enricher(
      gene = genes,
      TERM2GENE = term2gene,
      universe = universe,
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2
    ),
    error = function(e) NULL
  )

  if (is.null(result) || !nrow(as.data.frame(result))) return(data.frame())

  as.data.frame(result) %>%
    left_join(term2name, by = c("ID" = "GO_ID")) %>%
    mutate(Description = ifelse(!is.na(GO_Name), GO_Name, Description)) %>%
    select(-GO_Name) %>%
    arrange(p.adjust) %>%
    mutate(GeneRatio = as.character(GeneRatio), BgRatio = as.character(BgRatio))
}

run_kegg_enrichment <- function(genes, annotation, universe) {
  kegg_genes <- na.omit(unique(annotation$Old_Locus_Tag[match(genes, annotation$Locus_Tag)]))
  kegg_universe <- na.omit(unique(annotation$Old_Locus_Tag[match(universe, annotation$Locus_Tag)]))

  if (!length(kegg_genes) || !length(kegg_universe)) return(data.frame())

  result <- tryCatch(
    enrichKEGG(
      gene = kegg_genes,
      organism = kegg_organism,
      universe = kegg_universe,
      keyType = "kegg",
      pvalueCutoff = 0.05
    ),
    error = function(e) NULL
  )

  if (is.null(result) || !nrow(as.data.frame(result))) return(data.frame())

  as.data.frame(result) %>%
    arrange(p.adjust) %>%
    mutate(GeneRatio = as.character(GeneRatio), BgRatio = as.character(BgRatio))
}

make_summary_rows <- function(enrichment_df) {
  if (is.null(enrichment_df) || !nrow(enrichment_df)) return(data.frame())

  gene_total <- ratio_denominator(enrichment_df$GeneRatio)
  background_total <- ratio_denominator(enrichment_df$BgRatio)
  background_hits <- ratio_numerator(enrichment_df$BgRatio)

  data.frame(
    Cluster = enrichment_df$Cluster,
    Pathway = enrichment_df$Description,
    ProteinNumber = enrichment_df$Count,
    RichFactor = enrichment_df$Count / background_hits,
    FoldEnrichment = (enrichment_df$Count / gene_total) / (background_hits / background_total),
    padj = enrichment_df$p.adjust,
    stringsAsFactors = FALSE
  )
}

run_cluster_enrichment <- function(cluster_id, cluster_members, annotation, go_sets, universe) {
  all_members <- cluster_members$Locus_Tag[cluster_members$Cluster == cluster_id]
  core_members <- cluster_members$Locus_Tag[
    cluster_members$Cluster == cluster_id & cluster_members$Membership >= alpha_core
  ]

  sections <- list(
    KEGG = run_kegg_enrichment(all_members, annotation, universe),
    `GO BP - all members` = run_go_enrichment(all_members, go_sets$bp$term2gene, go_sets$bp$term2name, universe),
    `GO BP - membership >= 0.7` = run_go_enrichment(core_members, go_sets$bp$term2gene, go_sets$bp$term2name, universe),
    `GO Molecular Function - all members` = run_go_enrichment(all_members, go_sets$mf$term2gene, go_sets$mf$term2name, universe),
    `GO Molecular Function - membership >= 0.7` = run_go_enrichment(core_members, go_sets$mf$term2gene, go_sets$mf$term2name, universe)
  )

  bind_rows(lapply(names(sections), function(section_name) {
    section <- sections[[section_name]]
    if (!nrow(section)) return(data.frame())
    section$Cluster <- cluster_id
    section$Section <- section_name
    section
  }))
}


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

proteomics_results <- file.path(paths$output, "proteomics_mfuzz_results.xlsx")
if (!file.exists(proteomics_results)) {
  stop("Missing proteomics results.")
}

cluster_members <- openxlsx::read.xlsx(proteomics_results, sheet = "cluster_membership") %>%
  mutate(
    Cluster = as.integer(Cluster),
    Membership = as.numeric(Membership)
  )

annotation <- openxlsx::read.xlsx(paths$annotation) %>%
  mutate(Locus_Tag = as.character(Locus_Tag))

universe <- cluster_members$Locus_Tag

go_bp <- extract_go_terms(
  annotation %>% filter(Locus_Tag %in% universe),
  "GO-Biological_Process_ID",
  "GO-Biological_Process_Name"
)
go_mf <- extract_go_terms(
  annotation %>% filter(Locus_Tag %in% universe),
  "GO-Molecular_Function_ID",
  "GO-Molecular_Function_Name"
)

go_sets <- list(
  bp = list(
    term2gene = go_bp %>% select(GO_ID, Gene),
    term2name = go_bp %>% select(GO_ID, GO_Name) %>% distinct()
  ),
  mf = list(
    term2gene = go_mf %>% select(GO_ID, Gene),
    term2name = go_mf %>% select(GO_ID, GO_Name) %>% distinct()
  )
)

all_enrichments <- bind_rows(lapply(seq_len(n_clusters), function(cluster_id) {
  run_cluster_enrichment(cluster_id, cluster_members, annotation, go_sets, universe)
}))

summary <- all_enrichments %>%
  filter(Section %in% c("KEGG", "GO BP - all members", "GO Molecular Function - all members")) %>%
  make_summary_rows() %>%
  arrange(Cluster, padj)


# -----------------------------------------------------------------------------
# Write Results
# -----------------------------------------------------------------------------

out_file <- file.path(paths$output, "pathway_analysis_results.xlsx")

wb <- createWorkbook()
addWorksheet(wb, "summary")
writeData(wb, "summary", summary)
addWorksheet(wb, "all_enrichments")
writeData(wb, "all_enrichments", all_enrichments)
saveWorkbook(wb, out_file, overwrite = TRUE)

cat(sprintf(
  paste0(
    "Wrote %s\n",
    "Summary rows: %d\n",
    "Detailed enrichment rows: %d\n"
  ),
  normalizePath(out_file), nrow(summary), nrow(all_enrichments)
))