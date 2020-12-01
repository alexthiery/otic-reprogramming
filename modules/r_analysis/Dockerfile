FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run all R analysis for otic reprogramming project"

# Install apt packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    cabextract \
    libbz2-dev \
    liblzma-dev \
    build-essential=12.6 \
    libhdf5-dev=1.10.4+repack-10 \
    libboost-filesystem1.67-dev

# Install and unpack webfonts for plotting with Arial
RUN wget https://www.freedesktop.org/software/fontconfig/webfonts/webfonts.tar.gz && \
    tar -xzf webfonts.tar.gz && \
    cd msfonts/ && \
    cabextract *.exe && \
    cp *.ttf *.TTF /usr/share/fonts/

RUN R -e "install.packages('BiocManager')" && \
    R -e "BiocManager::install(version = '3.10')" && \
    R -e 'BiocManager::install("rtracklayer")' && \
    R -e 'BiocManager::install("ChIPpeakAnno")' && \
    R -e 'BiocManager::install("biomaRt")' && \
    R -e 'BiocManager::install("DESeq2")' && \
    R -e 'BiocManager::install("apeglm")' && \
    R -e "BiocManager::install('pcaMethods')" && \
    R -e "BiocManager::install('monocle')" && \
    R -e "devtools::install_github('velocyto-team/velocyto.R')" && \
    R -e "devtools::install_github('juliendelile/Antler', ref = 'Development2019')" && \
    R -e "devtools::install_github('hadley/multidplyr')" && \
    R -e 'install.packages(c("getopt", "ggsignif", "getopt", "plyr", "readr", "RColorBrewer", "VennDiagram", "pheatmap", "ggplot2", "openxlsx", "ggrepel", "ggseqlogo", "cowplot", "gridGraphics", "gprofiler2", "rstatix", "corrgram", "extrafont"))' && \
    # load fonts into R libs
    R -e 'library(extrafont); font_import(prompt = FALSE)'