#!/usr/bin/env Rscript

# Publication enrichment backend for the Network Pharmacology app.
# Mirrors the supplied EGCG enrichment notebooks while making the
# foreground/background choice explicit and machine-readable.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(
    paste(
      "Usage: Rscript run_enrichment.R",
      "<foreground.csv> <background_mode> <background.csv-or-NONE>",
      "<outdir> <fdr> <go_similarity> <reactome_similarity> <profile>"
    )
  )
}

foreground_file <- args[[1]]
background_mode <- args[[2]]
background_file <- args[[3]]
outdir <- args[[4]]
fdr_cutoff <- as.numeric(args[[5]])
go_similarity_cutoff <- as.numeric(args[[6]])
reactome_similarity_cutoff <- as.numeric(args[[7]])
analysis_profile <- args[[8]]

valid_background_modes <- c("package_default", "all_annotated", "custom")
if (!(background_mode %in% valid_background_modes)) {
  stop("Unsupported background mode: ", background_mode)
}
if (!(analysis_profile %in% c("publication", "legacy_notebook"))) {
  stop("Unsupported analysis profile: ", analysis_profile)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "AnnotationDbi",
  "GOSemSim",
  "jsonlite"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ReactomePA)
})

read_gene_symbols <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path)
  }
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(x) < 1) {
    stop("Input file contains no columns: ", path)
  }
  values <- trimws(as.character(x[[1]]))
  unique(values[!is.na(values) & values != ""])
}

map_symbols <- function(symbols) {
  if (length(symbols) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  }
  mapped <- suppressWarnings(
    clusterProfiler::bitr(
      symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  unique(mapped)
}

rank_enrichment <- function(df) {
  if (nrow(df) == 0) {
    df$Rank <- integer(0)
    return(df)
  }
  df <- df[order(df$p.adjust, -df$Count), , drop = FALSE]
  df$Rank <- seq_len(nrow(df))
  df[c("Rank", setdiff(colnames(df), "Rank"))]
}

filter_significant <- function(object) {
  df <- as.data.frame(object)
  df[
    !is.na(df$p.adjust) & df$p.adjust < fdr_cutoff,
    ,
    drop = FALSE
  ]
}

run_go <- function(ids, ontology, universe_ids = NULL) {
  pcut <- if (analysis_profile == "legacy_notebook") fdr_cutoff else 1
  qcut <- if (analysis_profile == "legacy_notebook") fdr_cutoff else 1

  if (is.null(universe_ids)) {
    clusterProfiler::enrichGO(
      gene = ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = pcut,
      qvalueCutoff = qcut,
      readable = TRUE
    )
  } else {
    clusterProfiler::enrichGO(
      gene = ids,
      universe = universe_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff = pcut,
      qvalueCutoff = qcut,
      readable = TRUE
    )
  }
}

run_reactome <- function(ids, universe_ids = NULL) {
  pcut <- if (analysis_profile == "legacy_notebook") fdr_cutoff else 1
  qcut <- if (analysis_profile == "legacy_notebook") fdr_cutoff else 1

  if (is.null(universe_ids)) {
    ReactomePA::enrichPathway(
      gene = ids,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = pcut,
      qvalueCutoff = qcut,
      readable = TRUE
    )
  } else {
    ReactomePA::enrichPathway(
      gene = ids,
      universe = universe_ids,
      organism = "human",
      pAdjustMethod = "BH",
      pvalueCutoff = pcut,
      qvalueCutoff = qcut,
      readable = TRUE
    )
  }
}

reduce_go <- function(enrichment_object, preselect_n = NULL) {
  df <- as.data.frame(enrichment_object)
  df <- df[
    !is.na(df$p.adjust) & df$p.adjust < fdr_cutoff,
    ,
    drop = FALSE
  ]

  if (nrow(df) == 0) {
    return(list(object = NULL, table = data.frame()))
  }

  if (!is.null(preselect_n) && nrow(df) > preselect_n) {
    df <- df[order(df$p.adjust, -df$Count), , drop = FALSE]
    df <- df[seq_len(preselect_n), , drop = FALSE]
  }

  ids <- df$ID
  subset_object <- enrichment_object
  subset_object@result <- enrichment_object@result[
    enrichment_object@result$ID %in% ids,
    ,
    drop = FALSE
  ]

  reduced_object <- clusterProfiler::simplify(
    subset_object,
    cutoff = go_similarity_cutoff,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang"
  )

  reduced_df <- as.data.frame(reduced_object)
  reduced_df <- reduced_df[
    !is.na(reduced_df$p.adjust) & reduced_df$p.adjust < fdr_cutoff,
    ,
    drop = FALSE
  ]

  if (nrow(reduced_df) > 0) {
    reduced_df <- reduced_df[
      order(reduced_df$p.adjust, -reduced_df$Count),
      ,
      drop = FALSE
    ]
    reduced_df$Rank <- seq_len(nrow(reduced_df))
    reduced_df <- reduced_df[
      c("Rank", setdiff(colnames(reduced_df), "Rank")),
      ,
      drop = FALSE
    ]
  }

  list(object = reduced_object, table = reduced_df)
}

reduce_reactome <- function(reactome_ranked) {
  reactome_sig <- reactome_ranked

  if (nrow(reactome_sig) == 0) {
    reactome_sig$RedundancyCluster <- integer(0)
    return(list(full = reactome_sig, reduced = reactome_sig))
  }

  if (nrow(reactome_sig) == 1) {
    reactome_sig$RedundancyCluster <- 1L
    reduced <- reactome_sig
    reduced$Rank <- 1L
    reduced <- reduced[c("Rank", setdiff(colnames(reduced), "Rank"))]
    return(list(full = reactome_sig, reduced = reduced))
  }

  pathway_genes <- strsplit(reactome_sig$geneID, "/", fixed = TRUE)
  names(pathway_genes) <- reactome_sig$ID

  jaccard <- function(x, y) {
    union_size <- length(union(x, y))
    if (union_size == 0) {
      return(0)
    }
    length(intersect(x, y)) / union_size
  }

  n_pathways <- length(pathway_genes)
  sim_matrix <- matrix(
    0,
    nrow = n_pathways,
    ncol = n_pathways,
    dimnames = list(names(pathway_genes), names(pathway_genes))
  )

  for (i in seq_len(n_pathways)) {
    for (j in i:n_pathways) {
      similarity <- jaccard(pathway_genes[[i]], pathway_genes[[j]])
      sim_matrix[i, j] <- similarity
      sim_matrix[j, i] <- similarity
    }
  }

  distance_matrix <- as.dist(1 - sim_matrix)
  reactome_hclust <- hclust(distance_matrix, method = "average")
  reactome_clusters <- cutree(
    reactome_hclust,
    h = 1 - reactome_similarity_cutoff
  )

  reactome_sig$RedundancyCluster <- reactome_clusters
  cluster_ids <- sort(unique(reactome_sig$RedundancyCluster))
  representative_rows <- integer(0)

  for (cluster_id in cluster_ids) {
    cluster_rows <- which(
      reactome_sig$RedundancyCluster == cluster_id
    )
    cluster_order <- order(
      reactome_sig$p.adjust[cluster_rows],
      -reactome_sig$Count[cluster_rows]
    )
    representative_rows <- c(
      representative_rows,
      cluster_rows[cluster_order[1]]
    )
  }

  reduced <- reactome_sig[
    representative_rows,
    ,
    drop = FALSE
  ]
  reduced <- reduced[
    order(reduced$p.adjust, -reduced$Count),
    ,
    drop = FALSE
  ]
  reduced$Rank <- seq_len(nrow(reduced))
  reduced <- reduced[c("Rank", setdiff(colnames(reduced), "Rank"))]

  list(full = reactome_sig, reduced = reduced)
}

foreground_symbols <- read_gene_symbols(foreground_file)
if (length(foreground_symbols) < 2) {
  stop("At least two foreground gene symbols are required.")
}

gene_map <- map_symbols(foreground_symbols)
foreground_ids <- unique(gene_map$ENTREZID)
missing_symbols <- setdiff(foreground_symbols, gene_map$SYMBOL)

if (length(foreground_ids) < 2) {
  stop("Fewer than two foreground genes mapped to Entrez IDs.")
}

human_universe <- NULL
background_label <- "Package/default enrichment universe"
background_mapping <- data.frame()

if (background_mode == "all_annotated") {
  human_universe <- unique(
    na.omit(
      AnnotationDbi::keys(
        org.Hs.eg.db,
        keytype = "ENTREZID"
      )
    )
  )
  background_label <- "All human Entrez IDs represented in org.Hs.eg.db"
}

if (background_mode == "custom") {
  if (background_file == "NONE" || !file.exists(background_file)) {
    stop("Custom background mode requires a background CSV.")
  }
  background_symbols <- read_gene_symbols(background_file)
  background_mapping <- map_symbols(background_symbols)
  human_universe <- unique(background_mapping$ENTREZID)
  background_label <- "User-supplied custom gene universe mapped through org.Hs.eg.db"
  if (length(human_universe) < 2) {
    stop("Custom background mapped to fewer than two Entrez IDs.")
  }
}

if (!is.null(human_universe)) {
  foreground_in_universe <- foreground_ids %in% human_universe
  if (any(!foreground_in_universe)) {
    outside <- foreground_ids[!foreground_in_universe]
    stop(
      "Foreground Entrez IDs outside the selected background: ",
      paste(outside, collapse = ", ")
    )
  }
}

go_bp <- run_go(foreground_ids, "BP", human_universe)
go_cc <- run_go(foreground_ids, "CC", human_universe)
go_mf <- run_go(foreground_ids, "MF", human_universe)
reactome <- run_reactome(foreground_ids, human_universe)

go_bp_ranked <- rank_enrichment(filter_significant(go_bp))
go_cc_ranked <- rank_enrichment(filter_significant(go_cc))
go_mf_ranked <- rank_enrichment(filter_significant(go_mf))
reactome_ranked <- rank_enrichment(filter_significant(reactome))

go_bp_reduction <- reduce_go(go_bp, preselect_n = 100)
go_cc_reduction <- reduce_go(go_cc)
go_mf_reduction <- reduce_go(go_mf)
reactome_reduction <- reduce_reactome(reactome_ranked)

go_bp_reduced <- go_bp_reduction$table
go_cc_reduced <- go_cc_reduction$table
go_mf_reduced <- go_mf_reduction$table
reactome_clustered <- reactome_reduction$full
reactome_reduced <- reactome_reduction$reduced

write.csv(gene_map, file.path(outdir, "symbol_to_entrez_mapping.csv"), row.names = FALSE)
write.csv(
  data.frame(GeneSymbol = missing_symbols),
  file.path(outdir, "unmapped_symbols.csv"),
  row.names = FALSE
)

if (!is.null(human_universe)) {
  write.csv(
    data.frame(ENTREZID = human_universe),
    file.path(outdir, "background_universe_entrez.csv"),
    row.names = FALSE
  )
}
if (nrow(background_mapping) > 0) {
  write.csv(
    background_mapping,
    file.path(outdir, "background_symbol_to_entrez_mapping.csv"),
    row.names = FALSE
  )
}

write.csv(go_bp_ranked, file.path(outdir, "go_bp_significant_ranked.csv"), row.names = FALSE)
write.csv(go_cc_ranked, file.path(outdir, "go_cc_significant_ranked.csv"), row.names = FALSE)
write.csv(go_mf_ranked, file.path(outdir, "go_mf_significant_ranked.csv"), row.names = FALSE)
write.csv(reactome_ranked, file.path(outdir, "reactome_significant_ranked.csv"), row.names = FALSE)

write.csv(go_bp_reduced, file.path(outdir, "go_bp_nonredundant.csv"), row.names = FALSE)
write.csv(go_cc_reduced, file.path(outdir, "go_cc_nonredundant.csv"), row.names = FALSE)
write.csv(go_mf_reduced, file.path(outdir, "go_mf_nonredundant.csv"), row.names = FALSE)
write.csv(reactome_clustered, file.path(outdir, "reactome_significant_clustered.csv"), row.names = FALSE)
write.csv(reactome_reduced, file.path(outdir, "reactome_nonredundant.csv"), row.names = FALSE)

save(
  go_bp,
  go_cc,
  go_mf,
  reactome,
  human_universe,
  foreground_ids,
  foreground_symbols,
  file = file.path(outdir, "enrichment_objects.RData")
)

package_versions <- data.frame(
  Package = c("R", required_packages),
  Version = c(
    R.version.string,
    vapply(
      required_packages,
      function(pkg) as.character(utils::packageVersion(pkg)),
      character(1)
    )
  )
)
write.csv(package_versions, file.path(outdir, "software_versions.csv"), row.names = FALSE)

summary <- list(
  status = "ok",
  organism = "Homo sapiens",
  analysis_profile = analysis_profile,
  background_mode = background_mode,
  background_definition = background_label,
  foreground_symbols = length(foreground_symbols),
  mapped_entrez_ids = length(foreground_ids),
  unmapped_symbols = as.list(missing_symbols),
  background_entrez_ids = if (is.null(human_universe)) NA else length(human_universe),
  fdr_cutoff = fdr_cutoff,
  multiple_testing = "Benjamini-Hochberg",
  go_similarity_measure = "Wang",
  go_similarity_cutoff = go_similarity_cutoff,
  go_bp_preselection_before_simplify = 100,
  reactome_similarity_measure = "Jaccard gene-set similarity",
  reactome_similarity_cutoff = reactome_similarity_cutoff,
  reactome_clustering = "average-linkage hierarchical clustering",
  significant_counts = list(
    go_bp = nrow(go_bp_ranked),
    go_cc = nrow(go_cc_ranked),
    go_mf = nrow(go_mf_ranked),
    reactome = nrow(reactome_ranked)
  ),
  nonredundant_counts = list(
    go_bp = nrow(go_bp_reduced),
    go_cc = nrow(go_cc_reduced),
    go_mf = nrow(go_mf_reduced),
    reactome = nrow(reactome_reduced)
  ),
  software_versions = as.list(
    stats::setNames(package_versions$Version, package_versions$Package)
  )
)

jsonlite::write_json(
  summary,
  file.path(outdir, "summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

methods <- c(
  "PUBLICATION ENRICHMENT ANALYSIS",
  "",
  paste0("Organism: Homo sapiens"),
  paste0("Background: ", background_label),
  paste0("Analysis profile: ", analysis_profile),
  paste0("Foreground symbols: ", length(foreground_symbols)),
  paste0("Mapped Entrez IDs: ", length(foreground_ids)),
  paste0("BH-adjusted P-value threshold: ", fdr_cutoff),
  paste0("GO redundancy reduction: Wang similarity cutoff = ", go_similarity_cutoff),
  "GO-BP preselection before simplify(): top 100 significant terms ranked by adjusted P, then Count.",
  paste0(
    "Reactome redundancy reduction: Jaccard cutoff = ",
    reactome_similarity_cutoff,
    "; average-linkage hierarchical clustering; representative = lowest adjusted P, then highest Count."
  ),
  "",
  "Raw significant tables and non-redundant interpretation tables are exported separately."
)
writeLines(methods, file.path(outdir, "METHODS.txt"))

cat(jsonlite::toJSON(summary, pretty = TRUE, auto_unbox = TRUE, na = "null"))
