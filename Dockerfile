# getting bas image
FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive
# ===================================
# Set maintainer
LABEL maintainer="Nima Nouri <nima.nouri@sanofi.com>" \
      description="PMCB SCB pipeline for processing scRNA-seq data"
LABEL version="1.0.0"
# ===================================
# installing python and all required packages
RUN apt-get update && \
        apt-get install -y \
        apt-utils \
        ed \
        less \
        locales \
        vim-tiny \
        wget \
        ca-certificates \
        apt-transport-https \
        curl \
        cmake \
        git \
        python3.9 \ 
        python3-pip \
        bzip2 \
        vim \
        nano \
        acl \
        software-properties-common \
        fonts-liberation \
        gcc \
        libssl-dev \
        libcurl4-openssl-dev \
        libhdf5-dev \
        libxml2-dev \
        libnlopt-dev \
        libicu-dev \
        libjpeg62 \
        libgeos-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        && \
        apt-get clean && \
        rm -rf /var/lib/apt/lists/*
# ===================================
# RUN pip3 install --upgrade pip
# ===================================
# Installing required PYTHON packages
# RUN pip3 install -r requirements.txt
RUN pip3 install \
  setuptools==65.6.3 \
  numpy==1.23.5 \
  pandas==1.5.2 \
  scipy==1.9.3 \
  matplotlib==3.6.2 \
  scikit-learn==1.2 \
  fa2==0.3.5 \
  h5py==3.7.0 \
  networkx==2.8.7 \
  scrublet
# ===================================
# Download and Install the latest version of pandoc from github (Update version as required)
ARG PANDOC_VERSION=2.19.2
RUN wget -q --show-progress --no-check-certificate https://github.com/jgm/pandoc/releases/download/${PANDOC_VERSION}/pandoc-${PANDOC_VERSION}-1-amd64.deb \
    && dpkg -i pandoc-${PANDOC_VERSION}-1-amd64.deb \
    && rm pandoc-${PANDOC_VERSION}-1-amd64.deb
# ===================================
#installing R and all required libraries
# update indices
RUN apt update -qq
# install two helper packages we need
RUN apt install apt-transport-https software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
# RUN apt-key adv --keyserver http://keyserver.ubuntu.com --recv-key 'E298A3A825C0D65DFD57CBB651716619E084DAB9'
RUN curl -sSL \
'http://keyserver.ubuntu.com/pks/lookup?op=get&search=0xE298A3A825C0D65DFD57CBB651716619E084DAB9' \
| apt-key add -
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN add-apt-repository 'deb http://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y \
    r-base \
    r-base-core \
		r-base-dev \
		r-recommended
# ===================================
# R packages
# install 'Rcpp' 1.0.9. The lastest 1.0.10 version (2023-01-22) is incompatible with 'reticulate'.
# see https://github.com/rstudio/reticulate/issues/1328
# RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Rcpp/Rcpp_1.0.9.tar.gz', \
#               repos=NULL, type='source')"
# Define R package installation parameters
ARG R_DEPS="c('devtools', 'optparse', 'dplyr', 'ggplot2', 'pheatmap', 'grid', \
            'gridExtra', 'cowplot', 'RColorBrewer', 'tidyr', 'versions', \
            'reticulate', 'purrr', 'kableExtra', 'RcppAnnoy', 'uwot', 'ggrepel', \
            'Seurat', 'igraph', 'gridtext', 'gplots', 'gtools', \
            'BiocManager', 'hdf5r', 'readxl', 'DT', 'data.tree', 'plotly', \
            'visNetwork', 'data.table', 'R.utils')"
RUN Rscript -e "install.packages(${R_DEPS}, clean=TRUE)"
# RE harmony error : https://github.com/immunogenomics/harmony/issues/166
# RUN Rscript -e "devtools::install_github('eddelbuettel/harmony', force=TRUE)"
# ===================================
ARG R_BIOC="c('biomaRt', 'Orthology.eg.db', 'org.Hs.eg.db', 'org.Mm.eg.db', \
              'limma','edgeR','ComplexHeatmap','Nebulosa','DelayedMatrixStats', \
              'slingshot')"
RUN Rscript -e "BiocManager::install(ask=FALSE)"
RUN Rscript -e "BiocManager::install(${R_BIOC})"
# ===================================
RUN Rscript -e "devtools::install_github('immunogenomics/harmony')"
RUN Rscript -e "devtools::install_github('immunogenomics/presto')"
RUN Rscript -e "devtools::install_github('Sanofi-Public/PMCB-Sargent')"
RUN Rscript -e "devtools::install_github('Sanofi-Public/PMCB-SignacX')"
# ===================================
RUN mkdir -p /opt/cellbridge_1.0.0
WORKDIR /opt/cellbridge_1.0.0
COPY ./* /opt/cellbridge_1.0.0/
RUN chmod +x Base.R
RUN ln Base.R /usr/local/bin/cellbridge
# ===================================
# RUN Rscript ./Mouse.R
# RUN Rscript ./MacacaFascicularis.R
# RUN Rscript ./HgncEnsembl.R
# ===================================
# The VOLUME instruction and the -v option to docker run tell docker to store 
# files in a directory on the host instead of in the container’s file system.
VOLUME /data
WORKDIR /data
# ===================================
# CMD ["/bin/bash"]