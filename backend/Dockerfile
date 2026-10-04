FROM bioconductor/bioconductor_docker:RELEASE_3_23

USER root

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-venv python3-pip curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY backend/requirements.txt /app/backend/requirements.txt
RUN python3 -m venv /opt/backend-venv \
    && /opt/backend-venv/bin/pip install --no-cache-dir -r /app/backend/requirements.txt

COPY r/enrichment_publication.R /app/r/enrichment_publication.R

RUN Rscript -e 'if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos="https://cloud.r-project.org"); BiocManager::install(c("clusterProfiler","ReactomePA","AnnotationDbi","GOSemSim","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db"), version="3.23", ask=FALSE, update=FALSE)' \
    && Rscript -e 'stopifnot(as.character(packageVersion("clusterProfiler"))=="4.20.0"); stopifnot(as.character(packageVersion("ReactomePA"))=="1.56.0"); stopifnot(as.character(packageVersion("org.Hs.eg.db"))=="3.23.1"); stopifnot(as.character(packageVersion("org.Mm.eg.db"))=="3.23.0"); stopifnot(as.character(packageVersion("org.Rn.eg.db"))=="3.23.0")'

COPY backend/app.py /app/backend/app.py

ENV PATH="/opt/backend-venv/bin:${PATH}"
ENV PYTHONUNBUFFERED=1
ENV ANALYSIS_TIMEOUT_SECONDS=1200

EXPOSE 10000

CMD ["sh", "-c", "uvicorn backend.app:app --host 0.0.0.0 --port ${PORT:-10000} --workers 1"]
