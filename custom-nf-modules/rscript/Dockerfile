FROM rocker/tidyverse:3.6.3

LABEL authors="alex.thiery@crick.ac.uk" \
      description="Docker image containing all requirements to run  analysis"

ARG WHEN

# Install apt packages
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential
    # cmake \
    # libbamtools-dev \
    # libboost-dev \
    # libboost-iostreams-dev \
    # libboost-log-dev \
    # libboost-system-dev \
    # libboost-test-dev \
    # libxml2-dev \
    # libz-dev \
    # curl \
    # libhdf5-cpp-100 \ 
    # libarmadillo7 \
    # libarmadillo-dev

#  git=1:2.20.1-2+deb10u3 \
#  apt-utils=1.8.2 \
#  unzip=6.0-23+deb10u1 \
#  procps=2:3.3.15-2 \
#  build-essential=12.6 \
#  zlib1g-dev=1:1.2.11.dfsg-1 \ 
#  libxt-dev \
#  libmagick++-dev \
#  python-pip

# RUN pip install pandas

# RUN R -e "devtools::install_github('velocyto-team/velocyto.R')"

# RUN R -e "devtools::install_github('juliendelile/Antler', dependencies = TRUE)"
# RUN R -e "install.packages(c('Rtsne', 'gProfileR', 'tidyverse', 'stringr', 'Seurat', 'BiocManager'))"
# RUN R -e "BiocManager::install('monocle')"
# RUN R -e "devtools::install_github('hadley/multidplyr')"
# RUN R -e "BiocManager::install('pcaMethods')"


# #' Extra packages are also required

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# # BiocManager::install("monocle")
# # devtools::install_github('hadley/multidplyr')

# # velocyto.R only runs on UNIX
# # BiocManager::install("pcaMethods")
# # devtools::install_github("velocyto-team/velocyto.R")


# RUN R -e "devtools::install_version('Seurat', version = '3.1.5', dependencies= NA)"
# RUN R -e "devtools::install_version('future', version = '1.17.0', dependencies= T)"
# RUN R -e "devtools::install_version('cowplot', version = '1.0.0', dependencies= T)"
# RUN R -e "devtools::install_version('clustree', version = '0.4.2', dependencies= NA)"
# RUN R -e "devtools::install_version('gridExtra', version = '2.3', dependencies=T)"
# RUN R -e "devtools::install_version('getopt', version = '1.20.3', dependencies=T)"
# RUN R -e "devtools::install_version('pheatmap', version = '1.0.12', dependencies=T)"
# RUN R -e "BiocManager::install('limma')"