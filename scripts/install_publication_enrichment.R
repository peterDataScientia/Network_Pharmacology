#!/usr/bin/env Rscript

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cran <- c("jsonlite")
for (pkg in cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

bioc <- c(
  "clusterProfiler",
  "ReactomePA",
  "AnnotationDbi",
  "GOSemSim",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Rn.eg.db"
)

missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

cat("\nPublication-enrichment R environment\n")
cat("R:", R.version.string, "\n")
for (pkg in c(cran, bioc)) {
  cat(pkg, ":", as.character(utils::packageVersion(pkg)), "\n")
}
