FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run downstream ChIP/ATAC/bulk-RNAseq"

# Install apt packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libjpeg-dev\
    libpng-dev \
    git=1:2.20.1-2+deb10u3 \
    apt-utils=1.8.2.1 \
    unzip=6.0-23+deb10u1 \
    procps=2:3.3.15-2 \
    build-essential=12.6 \
    libhdf5-dev=1.10.4+repack-10 \
    libboost-filesystem1.67-dev \
    zlib1g-dev=1:1.2.11.dfsg-1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.10')"

RUN R -e 'BiocManager::install("ChIPpeakAnno")'
RUN R -e 'BiocManager::install("biomaRt")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("apeglm")'

RUN R -e 'install.packages("getopt")'
RUN R -e 'install.packages("plyr")'
RUN R -e 'install.packages("readr")'
RUN R -e 'install.packages("RColorBrewer")'
RUN R -e 'install.packages("VennDiagram")'
RUN R -e 'install.packages("pheatmap")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("openxlsx")'
RUN R -e 'install.packages("ggrepel")'