FROM rocker/r-ver:4.6.1

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
    libfreetype6-dev libharfbuzz-dev libfribidi-dev \
    libpng-dev libtiff-dev libjpeg-dev libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages("BiocManager"); BiocManager::install(version="3.23", ask=FALSE); BiocManager::install(c("clusterProfiler","ReactomePA","AnnotationDbi","GOSemSim","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","jsonlite"), ask=FALSE, update=FALSE); stopifnot(as.character(packageVersion("clusterProfiler"))=="4.20.0"); stopifnot(as.character(packageVersion("ReactomePA"))=="1.56.0"); stopifnot(as.character(packageVersion("AnnotationDbi"))=="1.74.0"); stopifnot(as.character(packageVersion("GOSemSim"))=="2.38.3"); stopifnot(as.character(packageVersion("org.Hs.eg.db"))=="3.23.1"); stopifnot(as.character(packageVersion("org.Mm.eg.db"))=="3.23.0"); stopifnot(as.character(packageVersion("org.Rn.eg.db"))=="3.23.0")'

COPY r/enrichment_publication.R /app/r/enrichment_publication.R
COPY r/build_semantic_cache.R /app/r/build_semantic_cache.R

RUN mkdir -p /app/semantic_cache \
    && Rscript /app/r/build_semantic_cache.R /app/semantic_cache \
    && rm -rf /tmp/*

ENV PUBLICATION_SEMDATA_DIR=/app/semantic_cache
ENV PUBLICATION_USE_SEMDATA_CACHE=1

CMD ["R"]
