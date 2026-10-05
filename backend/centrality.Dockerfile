FROM rocker/r-ver:4.6.1

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /analysis

RUN apt-get update && apt-get install -y --no-install-recommends \
    libglpk-dev libxml2-dev libcurl4-openssl-dev libssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org")); install.packages(c("remotes","jsonlite")); remotes::install_version("igraph", version="2.3.4", repos="https://cloud.r-project.org", upgrade="never")'

COPY r/network_centrality.R /analysis/network_centrality.R

CMD ["R"]
