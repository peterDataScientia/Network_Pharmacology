#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: enrichment_publication.R <targets.txt> <outdir> <organism> <background_mode> <background_file_or_NONE>")
}

target_file <- args[[1]]
outdir <- args[[2]]
organism <- tolower(args[[3]])
background_mode <- tolower(args[[4]])
background_file <- args[[5]]

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
write.csv(data.frame(SYMBOL = unmapped_symbols), file.path(outdir, "unmapped_foreground.csv"), row.names = FALSE)

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
  write.csv(data.frame(ENTREZID = background_ids),
            file.path(outdir, "background_universe.csv"),
            row.names = FALSE)
}

safe_df <- function(x) {
  if (is.null(x)) return(data.frame())
  as.data.frame(x)
}

run_go <- function(ontology) {
  clusterProfiler::enrichGO(
    gene = foreground_ids,
    universe = background_ids,
    OrgDb = OrgDb,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
}

filter_sig <- function(df) {
  if (nrow(df) == 0 || !"p.adjust" %in% names(df)) return(df[0, , drop = FALSE])
  out <- df[!is.na(df$p.adjust) & df$p.adjust < 0.05, , drop = FALSE]
  if (nrow(out) > 0) {
    out <- out[order(out$p.adjust, -out$Count), , drop = FALSE]
  }
  out
}

reduce_go <- function(obj, preselect_n = NULL) {
  df <- filter_sig(safe_df(obj))
  if (nrow(df) == 0) return(df)

  # Match the reference notebook exactly: semantic reduction is performed
  # only on the statistically significant GO terms. GO-BP is additionally
  # restricted to the strongest 100 significant terms before simplify().
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

  reduced <- clusterProfiler::simplify(
    subset_obj,
    cutoff = 0.70,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang"
  )
  out <- filter_sig(safe_df(reduced))
  if (nrow(out) > 0) out <- out[order(out$p.adjust, -out$Count), , drop = FALSE]
  out
}

go_bp <- run_go("BP")
go_cc <- run_go("CC")
go_mf <- run_go("MF")

go_bp_raw <- filter_sig(safe_df(go_bp))
go_cc_raw <- filter_sig(safe_df(go_cc))
go_mf_raw <- filter_sig(safe_df(go_mf))

# Mirrors the uploaded EGCG notebook: preselect the strongest 100 BP terms
# before semantic simplification to control computational burden.
go_bp_reduced <- reduce_go(go_bp, preselect_n = 100)
go_cc_reduced <- reduce_go(go_cc)
go_mf_reduced <- reduce_go(go_mf)

reactome <- ReactomePA::enrichPathway(
  gene = foreground_ids,
  universe = background_ids,
  organism = organism,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  qvalueCutoff = 1,
  readable = TRUE
)
reactome_raw <- filter_sig(safe_df(reactome))

reduce_reactome <- function(df) {
  if (nrow(df) == 0) {
    df$RedundancyCluster <- integer(0)
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
  sim <- matrix(0, nrow = n, ncol = n,
                dimnames = list(names(pathway_genes), names(pathway_genes)))
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
    clusters <- cutree(hc, h = 0.30)
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

write.csv(go_bp_raw, file.path(outdir, "go_bp_raw_significant.csv"), row.names = FALSE)
write.csv(go_cc_raw, file.path(outdir, "go_cc_raw_significant.csv"), row.names = FALSE)
write.csv(go_mf_raw, file.path(outdir, "go_mf_raw_significant.csv"), row.names = FALSE)
write.csv(go_bp_reduced, file.path(outdir, "go_bp_reduced.csv"), row.names = FALSE)
write.csv(go_cc_reduced, file.path(outdir, "go_cc_reduced.csv"), row.names = FALSE)
write.csv(go_mf_reduced, file.path(outdir, "go_mf_reduced.csv"), row.names = FALSE)
write.csv(reactome_raw, file.path(outdir, "reactome_raw_significant.csv"), row.names = FALSE)
write.csv(reactome_reduced, file.path(outdir, "reactome_reduced.csv"), row.names = FALSE)

pkg_version <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(pkg))
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
  background_entrez_count = if (is.null(background_ids)) NA_integer_ else length(background_ids),
  significance = "BH-adjusted P < 0.05",
  go_redundancy = "Wang semantic similarity; cutoff 0.70; representative=min p.adjust",
  reactome_redundancy = "Jaccard gene-set similarity; cutoff 0.70; average-linkage h=0.30; representative=min p.adjust then max Count",
  counts = list(
    go_bp_raw = nrow(go_bp_raw),
    go_bp_reduced = nrow(go_bp_reduced),
    go_cc_raw = nrow(go_cc_raw),
    go_cc_reduced = nrow(go_cc_reduced),
    go_mf_raw = nrow(go_mf_raw),
    go_mf_reduced = nrow(go_mf_reduced),
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
