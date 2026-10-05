# Network Pharmacology & Toxicology Analyzer

A free, browser-based Streamlit workflow for researchers who already have a list of common targets and need to move from those targets to a transparent, reproducible network-analysis package.

## Current workflow

1. Paste or upload target identifiers.
2. Map identifiers with STRING for the selected organism.
3. Retrieve a functional or physical STRING PPI network at a user-selected confidence threshold, with no added neighbor nodes.
4. Display and export STRING's native high-resolution network figure (PNG and SVG).
5. Calculate four unweighted topology metrics on the connected filtered PPI: Degree, Betweenness, Closeness and Eigenvector centrality.
6. Rank the Top N targets for each metric and designate targets present in all four Top-N lists as **4/4 consensus hubs**.
7. Report mapped targets with no retained PPI edge separately rather than silently discarding them.
8. Retrieve GO Biological Process, GO Molecular Function, GO Cellular Component, KEGG and Reactome enrichment.
9. Download TSV tables, publication-quality figures and a complete ZIP package with analysis settings.

## Hub-target methodology

The primary hub workflow follows a consensus-centrality strategy:

```text
STRING PPI at selected confidence
        ↓
connected filtered topology
        ↓
Degree
Betweenness
Closeness
Eigenvector
        ↓
Top N from each ranking
        ↓
4/4 intersection
        ↓
Consensus hub targets
```

STRING confidence is used as an **edge-retention confidence threshold**. It is not interpreted as biochemical interaction strength, binding affinity or effect magnitude. The primary centrality analysis is therefore unweighted after the confidence filter is applied.

The default STRING threshold is 0.90 (highest confidence), but users may select other supported thresholds. The default consensus cutoff is Top 10 per metric.

A 4/4 consensus hub is a network-topology candidate. It is not, by itself, evidence of causality, target engagement or therapeutic importance.

## Scientific notes

- STRING mapping is performed before network/enrichment queries.
- No extra neighbor nodes are added to multi-target PPI queries.
- Isolated mapped targets are retained in downloadable results but excluded from topology-based hub ranking.
- If the connected network contains no more nodes than the selected Top N, consensus ranking has little discriminatory value; the app warns the user.
- Enrichment results are filtered by STRING-reported FDR.
- Results depend on the target list, organism, STRING version, interaction threshold and Top-N cutoff.

## Run locally

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Deploy free on Streamlit Community Cloud

1. Sign in to Streamlit Community Cloud with GitHub.
2. Create a new app from this repository.
3. Select branch `main` and entrypoint `app.py`.
4. Deploy.

No API key is required for STRING. The app sends submitted identifiers to the public STRING API.

## Planned extensions

Potential later modules include PPI-enrichment statistics, optional transcriptomic validation, additional topology algorithms, module detection and richer provenance/report generation.


## Publication Enrichment

The main application includes a publication-oriented enrichment mode alongside the existing quick STRING enrichment. The validated R/Bioconductor workflow is invoked from the **Enrichment** tab and does not disable or replace the STRING workflow.

Architecture:

```text
Streamlit/Python
      ↓
temporary foreground/background files
      ↓
Rscript r/enrichment_publication.R
      ↓
AnnotationDbi + organism OrgDb
clusterProfiler::enrichGO
ReactomePA::enrichPathway
      ↓
BH-adjusted P < 0.05
      ↓
raw significant tables
      +
GO Wang semantic reduction (0.70)
      +
Reactome Jaccard reduction (0.70)
      ↓
CSV + JSON + mirrored publication figure
```

The production UI is integrated into `app.py`; `publication_enrichment_dev.py` is retained as an isolated development harness.

### Background assistant

The prototype asks how the foreground was generated and supports three explicit enrichment universes:

- **Custom study background** — preferred when an experimental/panel-specific universe is known.
- **All annotated genes for organism** — explicit OrgDb Entrez universe; useful as a transparent general reference for database/prediction target lists.
- **Package/default background** — retained for reproduction of older analyses, not presented as a universal best choice.

### Reproducibility

The R engine records R, clusterProfiler, ReactomePA, AnnotationDbi, GOSemSim and organism annotation-package versions in `summary.json`.

Install the R environment locally with:

```bash
Rscript scripts/install_publication_enrichment.R
```

Run the EGCG–RISI validation with:

```bash
python validation/run_egcg_validation.py
```

The validation intentionally runs both the older default-background design and the newer explicit all-annotated-human background design as separate analyses.

A GitHub Actions workflow, `.github/workflows/validate-publication-enrichment.yml`, performs Python syntax checks, installs the Bioconductor environment and runs the EGCG validation. The EGCG reference workflow reproduced all eight raw/reduced term counts exactly before production integration.


## Publication Enrichment compute deployment

The heavy publication-grade GO/Reactome workflow runs outside the Streamlit
Community Cloud process. The production executor is **GitHub Actions** on this
public repository.

### Primary executor: GitHub Actions

Flow:

```text
Streamlit
   ↓ workflow_dispatch
GitHub Actions temporary Ubuntu runner
   ↓
validated R/Bioconductor container
   ↓
r/enrichment_publication.R
   ↓
1-day GitHub Actions result artifact
   ↓
Streamlit downloads and renders the result
```

The workflow is `.github/workflows/publication-enrichment-job.yml`. A
prebuilt validated container is maintained by
`.github/workflows/build-publication-image.yml` in GitHub Container Registry.
If the prebuilt image is temporarily unavailable, the job can build the same
`backend/Dockerfile` directly from the repository.

Every job validates the scientific environment before its result is accepted:

- R 4.6.1 / Bioconductor 3.23
- clusterProfiler 4.20.0
- ReactomePA 1.56.0
- AnnotationDbi 1.74.0
- GOSemSim 2.38.3
- organism-specific OrgDb version

For the exact 32-gene EGCG–RISI reference set with package/default background,
the workflow additionally requires exact historical parity:

```text
GO-BP      1432 -> 49
GO-CC        23 -> 17
GO-MF        54 -> 26
Reactome    276 -> 91
```

The app uses the original submitted gene symbols for publication enrichment and
maps SYMBOL → ENTREZ directly with AnnotationDbi. STRING preferred-name changes
therefore do not alter the enrichment foreground.

Configure Streamlit with a fine-grained GitHub token restricted to this
repository and **Actions: read and write**:

```toml
PUBLICATION_GITHUB_TOKEN = "<fine-grained GitHub token>"
```

Optional overrides are:

```toml
PUBLICATION_GITHUB_REPOSITORY = "peterDataScientia/Network_Pharmacology"
PUBLICATION_GITHUB_WORKFLOW = "publication-enrichment-job.yml"
PUBLICATION_GITHUB_REF = "main"
```

When `PUBLICATION_GITHUB_TOKEN` is configured, the deployed Streamlit app sends
Publication Enrichment jobs only to GitHub Actions. Identical successful
requests from the same code revision reuse their still-valid result artifact
instead of starting another job. A local R execution path is retained only for
local development and validation; there is no paid remote-compute fallback.

