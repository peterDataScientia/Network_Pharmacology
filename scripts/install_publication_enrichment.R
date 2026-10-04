#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
organism <- if (length(args) >= 1) tolower(args[[1]]) else "all"

user_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(normalizePath(user_lib, mustWork = FALSE), .libPaths()))
}

install_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org",
      lib = .libPaths()[1]
    )
  }
}

install_cran("BiocManager")
install_cran("jsonlite")

org_map <- c(
  human = "org.Hs.eg.db",
  mouse = "org.Mm.eg.db",
  rat = "org.Rn.eg.db"
)

if (organism == "all") {
  organism_packages <- unname(org_map)
} else if (organism %in% names(org_map)) {
  organism_packages <- unname(org_map[[organism]])
} else {
  stop("Organism must be one of: human, mouse, rat, all")
}

bioc <- unique(c(
  "clusterProfiler",
  "ReactomePA",
  "AnnotationDbi",
  "GOSemSim",
  organism_packages
))

missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  BiocManager::install(
    missing,
    ask = FALSE,
    update = FALSE,
    lib = .libPaths()[1]
  )
}

cat("\nPublication-enrichment R environment\n")
cat("R:", R.version.string, "\n")
cat("Library:", .libPaths()[1], "\n")
for (pkg in c("jsonlite", bioc)) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, ":", as.character(utils::packageVersion(pkg)), "\n")
  } else {
    cat(pkg, ": MISSING\n")
  }
}
