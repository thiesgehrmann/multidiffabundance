FROM mambaorg/micromamba

COPY mda /usr/bin/mda

RUN micromamba install \
      -y -n base \
      -c bioconda -c conda-forge \
      r-base \
      r-devtools \
      r-tidyverse \
      r-digest \
      r-lmerTest \
      r-reshape2  \
      r-vegan \
      r-GUniFrac \
      r-lme4 \
      r-matrix \
      bioconductor-phyloseq \
      bioconductor-ancombc && \
      eval "$(micromamba shell hook --shell bash)" && \
      micromamba activate base && \
    (echo "library(devtools); devtools::install_github('thiesgehrmann/multidiffabundance@devel', dependencies=FALSE)" | R --no-save) && \
    micromamba clean --all --yes
    
LABEL maintainer="Thies Gehrmann"

#ENTRYPOINT ["/usr/bin/mda"] # Do not provide entrypoint for now