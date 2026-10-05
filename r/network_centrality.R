#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(igraph)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: network_centrality.R <network.tsv> <outdir> <top_n>")
}

network_file <- args[[1]]
outdir <- args[[2]]
top_n <- as.integer(args[[3]])

if (is.na(top_n) || top_n < 1) {
  stop("top_n must be a positive integer")
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

edges <- read.delim(
  network_file,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

required_cols <- c("preferredName_A", "preferredName_B")
missing_cols <- setdiff(required_cols, colnames(edges))
if (length(missing_cols) > 0) {
  stop(paste("Missing required edge columns:", paste(missing_cols, collapse = ", ")))
}

edges <- edges[
  nzchar(trimws(edges$preferredName_A)) &
  nzchar(trimws(edges$preferredName_B)),
  ,
  drop = FALSE
]

if (nrow(edges) == 0) {
  write.csv(data.frame(), file.path(outdir, "centrality_all_genes.csv"), row.names = FALSE)
  write.csv(data.frame(), file.path(outdir, "centrality_consensus_rankings.csv"), row.names = FALSE)
  write.csv(data.frame(), file.path(outdir, "consensus_hubs_4of4.csv"), row.names = FALSE)
  write_json(
    list(
      status = "ok",
      nodes = 0,
      edges = 0,
      top_n_requested = top_n,
      top_n_effective = 0,
      hub_count = 0,
      metrics = c("Degree", "Betweenness", "Closeness", "Eigenvector"),
      engine = list(
        language = "R",
        R = R.version.string,
        package = "igraph",
        igraph = as.character(packageVersion("igraph"))
      )
    ),
    file.path(outdir, "centrality_summary.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  quit(status = 0)
}

edge_df <- data.frame(
  from = trimws(as.character(edges$preferredName_A)),
  to = trimws(as.character(edges$preferredName_B)),
  stringsAsFactors = FALSE
)

g <- graph_from_data_frame(edge_df, directed = FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

genes <- V(g)$name
n <- vcount(g)
m <- ecount(g)

degree_v <- degree(g, mode = "all", loops = FALSE, normalized = FALSE)
betweenness_v <- betweenness(
  g,
  v = V(g),
  directed = FALSE,
  weights = NA,
  normalized = TRUE
)
closeness_v <- closeness(
  g,
  vids = V(g),
  mode = "all",
  weights = NA,
  normalized = TRUE
)

eig <- tryCatch(
  eigen_centrality(
    g,
    directed = FALSE,
    scale = TRUE,
    weights = NA,
    options = list(maxiter = 5000)
  )$vector,
  error = function(e) {
    warning(conditionMessage(e))
    rep(NA_real_, n)
  }
)

# igraph returns a named vector; preserve the graph vertex order explicitly.
degree_v <- as.numeric(degree_v[genes])
betweenness_v <- as.numeric(betweenness_v[genes])
closeness_v <- as.numeric(closeness_v[genes])
eig <- as.numeric(eig[genes])

centrality <- data.frame(
  Gene = genes,
  Degree = degree_v,
  Betweenness = betweenness_v,
  Closeness = closeness_v,
  Eigenvector = eig,
  stringsAsFactors = FALSE
)

rank_desc_min <- function(x) {
  out <- rank(-x, ties.method = "min", na.last = "keep")
  if (anyNA(out)) {
    out[is.na(out)] <- nrow(centrality)
  }
  as.integer(out)
}

metrics <- c("Degree", "Betweenness", "Closeness", "Eigenvector")
for (metric in metrics) {
  centrality[[paste0(metric, " rank")]] <- rank_desc_min(centrality[[metric]])
}

centrality <- centrality[
  order(
    -centrality$Degree,
    -centrality$Betweenness,
    -centrality$Closeness,
    -centrality$Eigenvector,
    centrality$Gene,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]
row.names(centrality) <- NULL

effective_top_n <- min(top_n, nrow(centrality))
ranked <- centrality
membership_cols <- character()

for (metric in metrics) {
  ord <- order(-ranked[[metric]], ranked$Gene, na.last = TRUE)
  selected <- ranked$Gene[ord[seq_len(effective_top_n)]]
  col <- paste0("Top ", effective_top_n, " ", metric)
  ranked[[col]] <- ranked$Gene %in% selected
  membership_cols <- c(membership_cols, col)
}

ranked[["Consensus count"]] <- rowSums(ranked[, membership_cols, drop = FALSE])
ranked[["Consensus"]] <- paste0(ranked[["Consensus count"]], "/4")
ranked[["4/4 consensus hub"]] <- ranked[["Consensus count"]] == 4

rank_cols <- paste0(metrics, " rank")
ranked[["Mean rank"]] <- rowMeans(ranked[, rank_cols, drop = FALSE])

ranked <- ranked[
  order(
    -ranked[["Consensus count"]],
    ranked[["Mean rank"]],
    ranked[["Degree rank"]],
    ranked$Gene
  ),
  ,
  drop = FALSE
]
row.names(ranked) <- NULL

hubs <- ranked[ranked[["4/4 consensus hub"]], , drop = FALSE]
row.names(hubs) <- NULL

write.csv(
  centrality,
  file.path(outdir, "centrality_all_genes.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  ranked,
  file.path(outdir, "centrality_consensus_rankings.csv"),
  row.names = FALSE,
  na = ""
)
write.csv(
  hubs,
  file.path(outdir, "consensus_hubs_4of4.csv"),
  row.names = FALSE,
  na = ""
)

write_json(
  list(
    status = "ok",
    nodes = n,
    edges = m,
    top_n_requested = top_n,
    top_n_effective = effective_top_n,
    hub_count = nrow(hubs),
    metrics = metrics,
    engine = list(
      language = "R",
      R = R.version.string,
      package = "igraph",
      igraph = as.character(packageVersion("igraph"))
    )
  ),
  file.path(outdir, "centrality_summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
