# Publication Enrichment Prototype

This directory is an isolated development implementation of the planned
publication-grade enrichment workflow. It is **not used by the production
`main` branch** until validation passes.

## Scientific workflow

The `publication` profile mirrors the supplied all-human-background notebook:

1. clean unique human gene symbols;
2. map SYMBOL → ENTREZID with `clusterProfiler::bitr()` and `org.Hs.eg.db`;
3. record mapped and unmapped identifiers;
4. choose an explicit enrichment universe;
5. run GO-BP, GO-CC and GO-MF with `clusterProfiler::enrichGO()`;
6. run Reactome with `ReactomePA::enrichPathway()`;
7. use Benjamini-Hochberg multiple-testing correction;
8. retain significant results at adjusted P < 0.05;
9. rank by adjusted P ascending, then Count descending;
10. reduce GO redundancy with Wang semantic similarity (default cutoff 0.70);
11. preselect the strongest 100 significant GO-BP terms before `simplify()`,
    matching the supplied notebook;
12. reduce Reactome redundancy with Jaccard gene-set similarity, average-linkage
    hierarchical clustering, and a default 0.70 similarity threshold;
13. retain one Reactome representative per cluster using lowest adjusted P,
    then highest Count;
14. export raw significant and non-redundant tables separately.

The implementation also handles the edge case of exactly one significant
Reactome pathway without calling `hclust()`.

## Background modes

- `custom`: preferred when a defensible study-specific universe exists
  (e.g. genes passing RNA-seq detection/QC, measured proteome, assayed panel).
- `all_annotated`: all human Entrez IDs represented in `org.Hs.eg.db`.
  This is an explicit general reference for workflows that do not have an
  experimental detection universe.
- `package_default`: leaves `universe` unspecified. This mode is retained
  primarily to reproduce legacy analyses and should not be silently mixed
  with explicit-background results.

## Analysis profiles

- `publication`: sets enrichment p/q cutoffs to 1 and performs the final
  statistical decision explicitly from BH-adjusted P, matching the newer
  explicit-background notebook.
- `legacy_notebook`: reproduces the older executed notebook's p/q cutoff
  behavior and is used by the regression test.

## Install the R/Bioconductor backend

From the repository root:

```bash
Rscript publication_enrichment/install_packages.R
```

Required packages are:

- clusterProfiler
- org.Hs.eg.db
- ReactomePA
- AnnotationDbi
- GOSemSim
- jsonlite

## EGCG regression test

The included reference foreground contains the 32 EGCG–RISI shared targets
used in the supplied PPI/enrichment workflow.

Run:

```bash
python publication_enrichment/validate_egcg.py
```

The executed legacy notebook reported:

| Result | Significant | Non-redundant |
|---|---:|---:|
| GO-BP | 1432 | 49 |
| GO-CC | 23 | 17 |
| GO-MF | 54 | 26 |
| Reactome | 276 | 91 |

The validation script exits non-zero if these counts are not reproduced.
Because annotation databases change over time, a mismatch is not automatically
a software bug. The script prints the R/Bioconductor package versions and can
be run with `--allow-annotation-drift` for structural testing only.

The newer all-human-background notebook supplied for this project contains the
workflow but does not contain executed enrichment output for those cells, so
no all-human expected counts are hard-coded here. They must be generated and
reviewed before being treated as a reference result.

## Output contract

Each successful run writes:

```text
summary.json
METHODS.txt
software_versions.csv
symbol_to_entrez_mapping.csv
unmapped_symbols.csv
background_universe_entrez.csv            # explicit backgrounds only
background_symbol_to_entrez_mapping.csv   # custom background only
go_bp_significant_ranked.csv
go_cc_significant_ranked.csv
go_mf_significant_ranked.csv
reactome_significant_ranked.csv
go_bp_nonredundant.csv
go_cc_nonredundant.csv
go_mf_nonredundant.csv
reactome_significant_clustered.csv
reactome_nonredundant.csv
enrichment_objects.RData
```

Raw/significant and reduced/non-redundant results are intentionally separate.
The reduced tables are for concise biological interpretation; the significant
tables preserve the complete statistical result.

## Merge gate

Do not merge this feature to `main` until:

1. the R environment installs reproducibly;
2. the EGCG legacy regression test passes under the intended reference package
   versions, or any annotation-version change has been explicitly accepted;
3. the all-human-background run has been executed and inspected;
4. Streamlit handling of missing R/backend packages fails gracefully;
5. downloadable Methods/settings correctly record foreground, background,
   package versions, FDR, and redundancy parameters.
