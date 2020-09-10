FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run scRNAseq analysis in Antler-dev"

# Install apt packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential=12.6 \
    libhdf5-dev=1.10.4+repack-10 \
    libboost-filesystem1.67-dev

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.10')"
RUN R -e "BiocManager::install('pcaMethods')"
RUN R -e "devtools::install_github('velocyto-team/velocyto.R')"
RUN R -e "devtools::install_github('juliendelile/Antler_dev@994c1b5b', auth_token='538dd8b13c75e38ce275c12607877034a95ffb30')"
RUN R -e "devtools::install_github('hadley/multidplyr')"
RUN R -e "BiocManager::install('monocle')"

RUN R -e "devtools::install_version('Seurat', version = '3.1.5', dependencies= NA)"
RUN R -e "devtools::install_version('getopt', version = '1.20.3', dependencies=T)"

