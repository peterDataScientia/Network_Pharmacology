FROM bioconductor/bioconductor:3.23-R-4.6.1

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-pip python3-venv \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(c("clusterProfiler","ReactomePA","AnnotationDbi","GOSemSim","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","jsonlite"), ask=FALSE, update=FALSE)'

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir fastapi==0.118.0 uvicorn==0.37.0 pydantic==2.12.0

COPY r/enrichment_publication.R /app/r/enrichment_publication.R
COPY backend/app.py /app/backend/app.py

ENV PORT=10000
EXPOSE 10000

CMD ["sh","-c","uvicorn backend.app:app --host 0.0.0.0 --port ${PORT}"]
