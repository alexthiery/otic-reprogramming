FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run downstream ChIP/ATAC/bulk-RNAseq"

# Install apt packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential=12.6 \
    libhdf5-dev=1.10.4+repack-10 \
    libboost-filesystem1.67-dev


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