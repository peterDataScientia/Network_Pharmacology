from __future__ import annotations

import modal

APP_NAME = "network-pharmacology-publication"
SECRET_NAME = "network-pharmacology-publication-api"

app = modal.App(APP_NAME)

# Reuse the exact Docker environment already validated scientifically:
# R 4.6.1 / Bioconductor 3.23 / clusterProfiler 4.20.0 /
# ReactomePA 1.56.0 / org.Hs.eg.db 3.23.1.
image = modal.Image.from_dockerfile(
    "backend/Dockerfile",
    context_dir=".",
)

api_secret = modal.Secret.from_name(
    SECRET_NAME,
    required_keys=["PUBLICATION_API_KEY"],
)

@app.function(
    image=image,
    secrets=[api_secret],
    cpu=2.0,
    memory=4096,
    timeout=900,
    startup_timeout=900,
    min_containers=0,
    max_containers=2,
    scaledown_window=300,
)
@modal.asgi_app()
def web():
    from backend.app import app as fastapi_app
    return fastapi_app
