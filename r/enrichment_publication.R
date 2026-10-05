#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop(paste(
    "Usage: enrichment_publication.R <targets.txt> <outdir> <organism>",
    "<background_mode> <background_file_or_NONE> [analysis_settings_json_or_NONE]"
  ))
}

target_file <- args[[1]]
outdir <- args[[2]]
organism <- tolower(args[[3]])
background_mode <- tolower(args[[4]])
background_file <- args[[5]]
settings_file <- if (length(args) >= 6) args[[6]] else "NONE"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

required <- c("clusterProfiler", "ReactomePA", "AnnotationDbi", "GOSemSim", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(paste("Missing R packages:", paste(missing, collapse = ", ")))
}

org_map <- c(
  human = "org.Hs.eg.db",
  mouse = "org.Mm.eg.db",
  rat = "org.Rn.eg.db"
)
if (!organism %in% names(org_map)) {
  stop("Supported organisms: human, mouse, rat")
}
org_pkg <- unname(org_map[[organism]])
if (!requireNamespace(org_pkg, quietly = TRUE)) {
  stop(paste("Missing organism annotation package:", org_pkg))
}
OrgDb <- getExportedValue(org_pkg, org_pkg)

default_settings <- list(
  ontologies = c("BP", "CC", "MF"),
  include_reactome = TRUE,
  p_adjust_method = "BH",
  p_adjust_cutoff = 0.05,
  min_gs_size = 10L,
  max_gs_size = 500L,
  reduce_go = TRUE,
  go_similarity_cutoff = 0.70,
  go_bp_preselect_n = 100L,
  reduce_reactome = TRUE,
  reactome_jaccard_cutoff = 0.70
)

analysis_settings <- default_settings
if (
  settings_file != "NONE" &&
  nzchar(settings_file) &&
  file.exists(settings_file)
) {
  supplied <- jsonlite::fromJSON(settings_file, simplifyVector = TRUE)
  if (!is.list(supplied)) stop("Analysis settings JSON must contain an object.")
  for (name in names(supplied)) {
    analysis_settings[[name]] <- supplied[[name]]
  }
}

analysis_settings$ontologies <- unique(toupper(as.character(analysis_settings$ontologies)))
analysis_settings$ontologies <- intersect(
  analysis_settings$ontologies,
  c("BP", "CC", "MF")
)
analysis_settings$include_reactome <- isTRUE(analysis_settings$include_reactome)
analysis_settings$p_adjust_method <- as.character(analysis_settings$p_adjust_method)[1]
analysis_settings$p_adjust_cutoff <- as.numeric(analysis_settings$p_adjust_cutoff)[1]
analysis_settings$min_gs_size <- as.integer(analysis_settings$min_gs_size)[1]
analysis_settings$max_gs_size <- as.integer(analysis_settings$max_gs_size)[1]
analysis_settings$reduce_go <- isTRUE(analysis_settings$reduce_go)
analysis_settings$go_similarity_cutoff <- as.numeric(
  analysis_settings$go_similarity_cutoff
)[1]
analysis_settings$reduce_reactome <- isTRUE(analysis_settings$reduce_reactome)
analysis_settings$reactome_jaccard_cutoff <- as.numeric(
  analysis_settings$reactome_jaccard_cutoff
)[1]

bp_preselect <- analysis_settings$go_bp_preselect_n
if (is.null(bp_preselect) || length(bp_preselect) == 0 || is.na(bp_preselect[1])) {
  analysis_settings$go_bp_preselect_n <- NULL
} else {
  analysis_settings$go_bp_preselect_n <- as.integer(bp_preselect)[1]
}

if (
  length(analysis_settings$ontologies) == 0 &&
  !analysis_settings$include_reactome
) {
  stop("At least one enrichment database/ontology must be enabled.")
}
if (!analysis_settings$p_adjust_method %in% c("BH", "bonferroni", "holm", "BY")) {
  stop("Unsupported p_adjust_method.")
}
if (
  !is.finite(analysis_settings$p_adjust_cutoff) ||
  analysis_settings$p_adjust_cutoff <= 0 ||
  analysis_settings$p_adjust_cutoff > 1
) {
  stop("p_adjust_cutoff must be > 0 and <= 1.")
}
if (
  is.na(analysis_settings$min_gs_size) ||
  is.na(analysis_settings$max_gs_size) ||
  analysis_settings$min_gs_size < 1 ||
  analysis_settings$max_gs_size < analysis_settings$min_gs_size
) {
  stop("Invalid gene-set size range.")
}
if (
  !is.finite(analysis_settings$go_similarity_cutoff) ||
  analysis_settings$go_similarity_cutoff <= 0 ||
  analysis_settings$go_similarity_cutoff > 1
) {
  stop("Invalid GO semantic-similarity cutoff.")
}
if (
  !is.finite(analysis_settings$reactome_jaccard_cutoff) ||
  analysis_settings$reactome_jaccard_cutoff <= 0 ||
  analysis_settings$reactome_jaccard_cutoff > 1
) {
  stop("Invalid Reactome Jaccard cutoff.")
}
if (
  !is.null(analysis_settings$go_bp_preselect_n) &&
  analysis_settings$go_bp_preselect_n < 1
) {
  stop("go_bp_preselect_n must be positive or null.")
}

semantic_cache_dir <- Sys.getenv("PUBLICATION_SEMDATA_DIR", "/app/semantic_cache")
use_semantic_cache <- tolower(Sys.getenv("PUBLICATION_USE_SEMDATA_CACHE", "1")) %in%
  c("1", "true", "yes", "on")

semantic_cache_path <- function(ontology) {
  file.path(
    semantic_cache_dir,
    paste0(organism, "_", tolower(ontology), ".rds")
  )
}

get_semantic_data <- function(ontology) {
  if (!use_semantic_cache) return(NULL)
  path <- semantic_cache_path(ontology)
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

read_ids <- function(path) {
  x <- trimws(readLines(path, warn = FALSE))
  unique(x[nzchar(x)])
}

symbols <- read_ids(target_file)
if (length(symbols) == 0) stop("Foreground target list is empty.")

mapping <- AnnotationDbi::select(
  OrgDb,
  keys = symbols,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
)
mapping <- mapping[!is.na(mapping$ENTREZID) & nzchar(mapping$ENTREZID), , drop = FALSE]
mapping <- unique(mapping)
foreground_ids <- unique(as.character(mapping$ENTREZID))
mapped_symbols <- unique(as.character(mapping$SYMBOL))
unmapped_symbols <- setdiff(symbols, mapped_symbols)

if (length(foreground_ids) < 2) {
  stop("Fewer than two foreground genes mapped to Entrez IDs.")
}

write.csv(mapping, file.path(outdir, "mapping_foreground.csv"), row.names = FALSE)
write.csv(
  data.frame(SYMBOL = unmapped_symbols),
  file.path(outdir, "unmapped_foreground.csv"),
  row.names = FALSE
)

background_ids <- NULL
background_description <- "package/default universe"

if (background_mode == "annotated") {
  background_ids <- unique(na.omit(AnnotationDbi::keys(OrgDb, keytype = "ENTREZID")))
  background_description <- paste("all Entrez IDs represented in", org_pkg)
} else if (background_mode == "custom") {
  if (background_file == "NONE" || !file.exists(background_file)) {
    stop("Custom background mode requires a background file.")
  }
  bg_symbols <- read_ids(background_file)
  bg_map <- AnnotationDbi::select(
    OrgDb,
    keys = bg_symbols,
    keytype = "SYMBOL",
    columns = c("SYMBOL", "ENTREZID")
  )
  bg_map <- bg_map[!is.na(bg_map$ENTREZID) & nzchar(bg_map$ENTREZID), , drop = FALSE]
  bg_map <- unique(bg_map)
  background_ids <- unique(as.character(bg_map$ENTREZID))
  background_description <- "user-supplied custom background mapped to Entrez IDs"
  write.csv(bg_map, file.path(outdir, "mapping_background.csv"), row.names = FALSE)
} else if (background_mode != "default") {
  stop("background_mode must be one of: default, annotated, custom")
}

if (!is.null(background_ids)) {
  outside <- setdiff(foreground_ids, background_ids)
  if (length(outside) > 0) {
    stop(paste(
      "Foreground Entrez IDs outside selected background:",
      paste(outside, collapse = ", ")
    ))
  }
  write.csv(
    data.frame(ENTREZID = background_ids),
    file.path(outdir, "background_universe.csv"),
    row.names = FALSE
  )
}

safe_df <- function(x) {
  if (is.null(x)) return(data.frame())
  as.data.frame(x)
}

run_go <- function(ontology) {
  if (!ontology %in% analysis_settings$ontologies) return(NULL)
  clusterProfiler::enrichGO(
    gene = foreground_ids,
    universe = background_ids,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = analysis_settings$p_adjust_method,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = analysis_settings$min_gs_size,
    maxGSSize = analysis_settings$max_gs_size,
    readable = TRUE
  )
}

filter_sig <- function(df) {
  if (nrow(df) == 0 || !"p.adjust" %in% names(df)) {
    return(df[0, , drop = FALSE])
  }
  out <- df[
    !is.na(df$p.adjust) &
      df$p.adjust < analysis_settings$p_adjust_cutoff,
    ,
    drop = FALSE
  ]
  if (nrow(out) > 0) {
    out <- out[order(out$p.adjust, -out$Count), , drop = FALSE]
  }
  out
}

reduce_go <- function(obj, ontology, preselect_n = NULL) {
  df <- filter_sig(safe_df(obj))
  if (nrow(df) == 0) return(df)
  if (!analysis_settings$reduce_go) return(df)

  keep_ids <- df$ID
  if (!is.null(preselect_n) && length(keep_ids) > preselect_n) {
    keep_ids <- head(keep_ids, preselect_n)
  }

  subset_obj <- obj
  subset_obj@result <- obj@result[
    obj@result$ID %in% keep_ids,
    ,
    drop = FALSE
  ]

  semdata <- get_semantic_data(ontology)

  reduced <- clusterProfiler::simplify(
    subset_obj,
    cutoff = analysis_settings$go_similarity_cutoff,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang",
    semData = semdata
  )
  out <- filter_sig(safe_df(reduced))
  if (nrow(out) > 0) {
    out <- out[order(out$p.adjust, -out$Count), , drop = FALSE]
  }
  out
}

go_bp <- run_go("BP")
go_cc <- run_go("CC")
go_mf <- run_go("MF")

go_bp_all <- safe_df(go_bp)
go_cc_all <- safe_df(go_cc)
go_mf_all <- safe_df(go_mf)

go_bp_raw <- filter_sig(go_bp_all)
go_cc_raw <- filter_sig(go_cc_all)
go_mf_raw <- filter_sig(go_mf_all)

go_bp_reduced <- reduce_go(
  go_bp,
  "BP",
  preselect_n = analysis_settings$go_bp_preselect_n
)
go_cc_reduced <- reduce_go(go_cc, "CC")
go_mf_reduced <- reduce_go(go_mf, "MF")

if (analysis_settings$include_reactome) {
  reactome <- ReactomePA::enrichPathway(
    gene = foreground_ids,
    universe = background_ids,
    organism = organism,
    pAdjustMethod = analysis_settings$p_adjust_method,
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = analysis_settings$min_gs_size,
    maxGSSize = analysis_settings$max_gs_size,
    readable = TRUE
  )
} else {
  reactome <- NULL
}

reactome_all <- safe_df(reactome)
reactome_raw <- filter_sig(reactome_all)

reduce_reactome <- function(df) {
  if (nrow(df) == 0) {
    df$RedundancyCluster <- integer(0)
    return(df)
  }
  if (!analysis_settings$reduce_reactome) {
    df$RedundancyCluster <- seq_len(nrow(df))
    return(df)
  }

  pathway_genes <- strsplit(as.character(df$geneID), "/", fixed = TRUE)
  names(pathway_genes) <- df$ID

  jaccard <- function(x, y) {
    denom <- length(union(x, y))
    if (denom == 0) return(0)
    length(intersect(x, y)) / denom
  }

  n <- length(pathway_genes)
  sim <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(names(pathway_genes), names(pathway_genes))
  )
  for (i in seq_len(n)) {
    for (j in i:n) {
      value <- jaccard(pathway_genes[[i]], pathway_genes[[j]])
      sim[i, j] <- value
      sim[j, i] <- value
    }
  }

  if (n == 1) {
    clusters <- 1L
  } else {
    hc <- hclust(as.dist(1 - sim), method = "average")
    # Preserve exact boundary semantics from the validated reference workflow.
    # Direct subtraction can produce e.g. 1 - 0.70 = 0.30000000000000004,
    # which can merge a cluster exactly at distance 0.30 that the historical
    # literal h = 0.30 kept separate.
    cluster_height <- round(
      1 - analysis_settings$reactome_jaccard_cutoff,
      digits = 12
    )
    clusters <- cutree(
      hc,
      h = cluster_height
    )
  }

  df$RedundancyCluster <- unname(clusters)
  reps <- integer(0)
  for (cluster_id in sort(unique(df$RedundancyCluster))) {
    rows <- which(df$RedundancyCluster == cluster_id)
    ord <- order(df$p.adjust[rows], -df$Count[rows])
    reps <- c(reps, rows[ord[[1]]])
  }
  out <- df[reps, , drop = FALSE]
  out[order(out$p.adjust, -out$Count), , drop = FALSE]
}

reactome_reduced <- reduce_reactome(reactome_raw)

write.csv(go_bp_all, file.path(outdir, "go_bp_all_tested.csv"), row.names = FALSE)
write.csv(go_cc_all, file.path(outdir, "go_cc_all_tested.csv"), row.names = FALSE)
write.csv(go_mf_all, file.path(outdir, "go_mf_all_tested.csv"), row.names = FALSE)
write.csv(reactome_all, file.path(outdir, "reactome_all_tested.csv"), row.names = FALSE)

write.csv(go_bp_raw, file.path(outdir, "go_bp_raw_significant.csv"), row.names = FALSE)
write.csv(go_cc_raw, file.path(outdir, "go_cc_raw_significant.csv"), row.names = FALSE)
write.csv(go_mf_raw, file.path(outdir, "go_mf_raw_significant.csv"), row.names = FALSE)
write.csv(go_bp_reduced, file.path(outdir, "go_bp_reduced.csv"), row.names = FALSE)
write.csv(go_cc_reduced, file.path(outdir, "go_cc_reduced.csv"), row.names = FALSE)
write.csv(go_mf_reduced, file.path(outdir, "go_mf_reduced.csv"), row.names = FALSE)
write.csv(
  reactome_raw,
  file.path(outdir, "reactome_raw_significant.csv"),
  row.names = FALSE
)
write.csv(
  reactome_reduced,
  file.path(outdir, "reactome_reduced.csv"),
  row.names = FALSE
)

pkg_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
}

bp_preselect_summary <- if (is.null(analysis_settings$go_bp_preselect_n)) {
  "all"
} else {
  as.integer(analysis_settings$go_bp_preselect_n)
}

summary <- list(
  status = "ok",
  organism = organism,
  organism_db = org_pkg,
  foreground_symbols_submitted = length(symbols),
  foreground_entrez_mapped = length(foreground_ids),
  unmapped_symbols = unmapped_symbols,
  background_mode = background_mode,
  background_description = background_description,
  background_entrez_count = if (is.null(background_ids)) {
    NA_integer_
  } else {
    length(background_ids)
  },
  significance = paste0(
    analysis_settings$p_adjust_method,
    "-adjusted P < ",
    format(analysis_settings$p_adjust_cutoff, scientific = FALSE)
  ),
  analysis_settings = list(
    ontologies = analysis_settings$ontologies,
    include_reactome = analysis_settings$include_reactome,
    p_adjust_method = analysis_settings$p_adjust_method,
    p_adjust_cutoff = analysis_settings$p_adjust_cutoff,
    min_gs_size = analysis_settings$min_gs_size,
    max_gs_size = analysis_settings$max_gs_size,
    reduce_go = analysis_settings$reduce_go,
    go_similarity_cutoff = analysis_settings$go_similarity_cutoff,
    go_bp_preselect_n = bp_preselect_summary,
    reduce_reactome = analysis_settings$reduce_reactome,
    reactome_jaccard_cutoff = analysis_settings$reactome_jaccard_cutoff
  ),
  go_redundancy = if (analysis_settings$reduce_go) {
    paste0(
      "Wang semantic similarity; cutoff ",
      analysis_settings$go_similarity_cutoff,
      "; representative=min p.adjust"
    )
  } else {
    "disabled"
  },
  go_semantic_cache = if (use_semantic_cache) {
    "precomputed when available"
  } else {
    "disabled"
  },
  reactome_redundancy = if (analysis_settings$reduce_reactome) {
    paste0(
      "Jaccard gene-set similarity; cutoff ",
      analysis_settings$reactome_jaccard_cutoff,
      "; average-linkage; representative=min p.adjust then max Count"
    )
  } else {
    "disabled"
  },
  counts = list(
    go_bp_tested = nrow(go_bp_all),
    go_bp_raw = nrow(go_bp_raw),
    go_bp_reduced = nrow(go_bp_reduced),
    go_cc_tested = nrow(go_cc_all),
    go_cc_raw = nrow(go_cc_raw),
    go_cc_reduced = nrow(go_cc_reduced),
    go_mf_tested = nrow(go_mf_all),
    go_mf_raw = nrow(go_mf_raw),
    go_mf_reduced = nrow(go_mf_reduced),
    reactome_tested = nrow(reactome_all),
    reactome_raw = nrow(reactome_raw),
    reactome_reduced = nrow(reactome_reduced)
  ),
  versions = list(
    R = R.version.string,
    clusterProfiler = pkg_version("clusterProfiler"),
    ReactomePA = pkg_version("ReactomePA"),
    AnnotationDbi = pkg_version("AnnotationDbi"),
    GOSemSim = pkg_version("GOSemSim"),
    organism_db = pkg_version(org_pkg)
  )
)

jsonlite::write_json(
  summary,
  file.path(outdir, "summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE,
  na = "null"
)

cat("Publication enrichment completed successfully\n")
