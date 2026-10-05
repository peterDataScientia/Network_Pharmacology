FROM bioconductor/bioconductor:3.23-R-4.6.1

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager"); BiocManager::install(c("clusterProfiler","ReactomePA","AnnotationDbi","GOSemSim","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","jsonlite"), ask=FALSE, update=FALSE); if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes"); remotes::install_version("igraph", version="2.3.4", repos="https://cloud.r-project.org", upgrade="never")'

COPY r/enrichment_publication.R /app/r/enrichment_publication.R
COPY r/network_centrality.R /app/r/network_centrality.R

CMD ["R"]
