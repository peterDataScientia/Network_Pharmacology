#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[[1]] else "/app/semantic_cache"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

organisms <- list(
  human = list(pkg = "org.Hs.eg.db"),
  mouse = list(pkg = "org.Mm.eg.db"),
  rat = list(pkg = "org.Rn.eg.db")
)

for (organism in names(organisms)) {
  pkg <- organisms[[organism]]$pkg
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste("Missing organism package:", pkg))
  }
  OrgDb <- getExportedValue(pkg, pkg)

  for (ontology in c("BP", "CC", "MF")) {
    message("Building semantic cache: ", organism, " ", ontology)
    semdata <- GOSemSim::godata(
      OrgDb = OrgDb,
      ont = ontology,
      computeIC = FALSE
    )
    saveRDS(
      semdata,
      file.path(outdir, paste0(organism, "_", tolower(ontology), ".rds")),
      compress = "xz"
    )
  }
}
