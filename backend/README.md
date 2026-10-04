# Publication Enrichment Backend

This service runs the publication-grade enrichment workflow outside Streamlit so users never compile Bioconductor packages during an interactive session.

## Scientific runtime

The Docker image is built from the official Bioconductor 3.23 container and the build fails unless these validated versions are present:

- R 4.6 / Bioconductor 3.23
- clusterProfiler 4.20.0
- ReactomePA 1.56.0
- org.Hs.eg.db 3.23.1
- org.Mm.eg.db 3.23.0
- org.Rn.eg.db 3.23.0

The same `r/enrichment_publication.R` script used in repository validation is copied into the image.

## API

- `GET /health` — runtime/package readiness and versions
- `POST /v1/enrichment` — GO/Reactome publication enrichment

When `BACKEND_API_TOKEN` is set, requests to the analysis endpoint must send:

```text
Authorization: Bearer <token>
```

Example request:

```json
{
  "targets": ["AKT1", "TP53", "STAT3"],
  "organism": "human",
  "background_mode": "annotated",
  "custom_background": []
}
```

Only one publication analysis is accepted at a time on a service instance. A concurrent request receives HTTP 429 instead of competing for memory.

## Render deployment

The repository root contains `render.yaml`. In Render, create a Blueprint from this repository and enter a private value for `BACKEND_API_TOKEN`.

Then add these root-level secrets to the Streamlit app:

```toml
PUBLICATION_BACKEND_URL = "https://<your-service>.onrender.com"
PUBLICATION_BACKEND_TOKEN = "<same-token>"
```

The Streamlit app checks `/health` before enabling publication analysis.

The free Render plan can be used for testing, but it has 512 MB RAM, 0.1 CPU, and spins down after inactivity. The exact workflow is CPU-heavy; a larger service will have much better latency. The backend architecture is independent of Render and can be deployed unchanged on another Docker host.
