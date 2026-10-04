options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

BiocManager::install(
  c(
    "clusterProfiler",
    "org.Hs.eg.db",
    "ReactomePA",
    "AnnotationDbi",
    "GOSemSim"
  ),
  ask = FALSE,
  update = FALSE
)

cat("Publication-enrichment packages installed.\n")
