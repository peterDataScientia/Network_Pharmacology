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

Potential later modules include PPI-enrichment statistics, R/Bioconductor enrichment with redundancy reduction, optional transcriptomic validation, additional topology algorithms, module detection and richer provenance/report generation.
