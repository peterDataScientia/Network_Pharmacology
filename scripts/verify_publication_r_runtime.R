#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
description_file <- if (length(args) >= 1) args[[1]] else "DESCRIPTION"

if (!file.exists(description_file)) {
  stop("Dependency manifest not found: ", description_file)
}

desc <- read.dcf(description_file)[1, ]

parse_dependencies <- function(field) {
  raw <- desc[[field]]
  if (is.null(raw) || is.na(raw) || !nzchar(raw)) return(data.frame())

  entries <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  entries <- entries[nzchar(entries)]

  rows <- lapply(entries, function(entry) {
    m <- regexec(
      "^([A-Za-z0-9.]+)(?:[[:space:]]*\\((==|>=|<=|>|<)[[:space:]]*([^)]+)\\))?$",
      entry
    )
    hit <- regmatches(entry, m)[[1]]
    if (length(hit) == 0) stop("Cannot parse dependency entry: ", entry)

    data.frame(
      package = hit[[2]],
      operator = if (length(hit) >= 3) hit[[3]] else "",
      version = if (length(hit) >= 4) trimws(hit[[4]]) else "",
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

satisfies <- function(got, op, expected) {
  if (!nzchar(op) || !nzchar(expected)) return(TRUE)
  cmp <- utils::compareVersion(got, expected)
  switch(
    op,
    "==" = cmp == 0,
    ">=" = cmp >= 0,
    "<=" = cmp <= 0,
    ">"  = cmp > 0,
    "<"  = cmp < 0,
    stop("Unsupported version operator: ", op)
  )
}

depends <- parse_dependencies("Depends")
imports <- parse_dependencies("Imports")

r_req <- depends[depends$package == "R", , drop = FALSE]
if (nrow(r_req) > 0) {
  got_r <- paste(R.version$major, R.version$minor, sep = ".")
  if (!satisfies(got_r, r_req$operator[[1]], r_req$version[[1]])) {
    stop(
      "R version mismatch: manifest requires R ",
      r_req$operator[[1]], " ", r_req$version[[1]],
      "; runtime has ", got_r
    )
  }
  cat("R\t", got_r, "\n", sep = "")
}

if (nrow(imports) == 0) stop("No Imports declared in DESCRIPTION.")

for (i in seq_len(nrow(imports))) {
  pkg <- imports$package[[i]]
  op <- imports$operator[[i]]
  expected <- imports$version[[i]]

  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg)
  }

  got <- as.character(utils::packageVersion(pkg))
  if (!satisfies(got, op, expected)) {
    stop(
      pkg, " version mismatch: manifest requires ",
      if (nzchar(op)) paste(op, expected) else "an installed version",
      "; runtime has ", got
    )
  }

  cat(pkg, "\t", got, "\n", sep = "")
}

cat("DECLARED_R_RUNTIME_OK\n")
