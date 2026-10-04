options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = max(1L, min(2L, parallel::detectCores())))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}

# The supplied notebooks use the R 4.6 generation. Bioconductor 3.23 is the
# matching release and is pinned here for reproducibility.
BiocManager::install(version = "3.23", ask = FALSE, update = FALSE)

required_bioc <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "ReactomePA",
  "AnnotationDbi",
  "GOSemSim"
)

BiocManager::install(
  required_bioc,
  ask = FALSE,
  update = FALSE,
  dependencies = TRUE
)

required_all <- c(required_bioc, "jsonlite")
missing <- required_all[
  !vapply(required_all, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  cat(
    "ERROR: required publication-enrichment package(s) are still missing: ",
    paste(missing, collapse = ", "),
    "\n",
    sep = ""
  )
  quit(save = "no", status = 1)
}

cat("Publication-enrichment packages verified successfully.\n")
for (pkg in required_all) {
  cat(pkg, ": ", as.character(packageVersion(pkg)), "\n", sep = "")
}
