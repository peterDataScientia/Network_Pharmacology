# Network Pharmacology & Toxicology Analyzer

A free, browser-based Streamlit workflow for researchers who already have a list of common targets and need to move from those targets to a reproducible network-analysis package.

## MVP workflow

1. Paste or upload target identifiers.
2. Map identifiers with STRING for the selected organism.
3. Retrieve a functional or physical STRING PPI network at a user-selected confidence threshold.
4. Calculate local centrality metrics and rank candidate hub genes.
5. Retrieve GO Biological Process, GO Molecular Function, GO Cellular Component, KEGG and Reactome enrichment.
6. Download TSV tables, 600-dpi PNG figures, PDF/SVG vector outputs and a complete ZIP package with analysis settings.

## Scientific notes

- STRING mapping is performed before network/enrichment queries.
- No extra neighbor nodes are added to the submitted PPI network.
- For shortest-path centrality, inverse STRING confidence is used as edge distance; STRING confidence itself is used as edge strength for weighted metrics.
- Enrichment results are filtered by STRING-reported FDR.
- Hub rankings are topology-dependent and should not be treated as causal biological evidence.

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

## Current scope

The first version intentionally stays lightweight enough for free hosting. A later version can add R/Bioconductor workflows, alternative enrichment backends, Cytoscape-style hub algorithms, provenance manifests, and richer interactive network editing.
