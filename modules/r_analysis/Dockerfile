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
RUN wget https://www.freedesktop.org/software/fontconfig/webfonts/webfonts.tar.gz &&\
    tar -xzf webfonts.tar.gz && \
    cd msfonts/ && \
    cabextract *.exeD && \
    cp *.ttf *.TTF /usr/share/fonts/

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.10')"

RUN R -e 'BiocManager::install("rtracklayer")'
RUN R -e 'BiocManager::install("ChIPpeakAnno")'
RUN R -e 'BiocManager::install("biomaRt")'
RUN R -e 'BiocManager::install("DESeq2")'
RUN R -e 'BiocManager::install("apeglm")'
RUN R -e "BiocManager::install('pcaMethods')"
RUN R -e "BiocManager::install('monocle')"

RUN R -e "devtools::install_github('velocyto-team/velocyto.R')"
RUN R -e "devtools::install_github('juliendelile/Antler', ref = 'Development2019')"
RUN R -e "devtools::install_github('hadley/multidplyr')"

RUN R -e 'install.packages("getopt")'
RUN R -e 'install.packages("ggsignif")'
RUN R -e 'install.packages("getopt")'
RUN R -e 'install.packages("plyr")'
RUN R -e 'install.packages("readr")'
RUN R -e 'install.packages("RColorBrewer")'
RUN R -e 'install.packages("VennDiagram")'
RUN R -e 'install.packages("pheatmap")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("openxlsx")'
RUN R -e 'install.packages("ggrepel")'
RUN R -e 'install.packages("ggseqlogo")'
RUN R -e 'install.packages("cowplot")'
RUN R -e 'install.packages("gridGraphics")'
RUN R -e 'install.packages("gprofiler2")'
RUN R -e 'install.packages("rstatix")'
RUN R -e 'install.packages("corrgram")'