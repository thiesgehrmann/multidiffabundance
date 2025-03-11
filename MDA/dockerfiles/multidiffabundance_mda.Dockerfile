FROM mambaorg/micromamba

COPY mda /usr/bin/mda

RUN micromamba install \
      -y -n base \
      -c bioconda -c conda-forge \
      cmake \
      r-base \
      r-devtools \
      r-tidyverse \
      r-digest \
      r-reshape2  \
      r-vegan \
      r-GUniFrac  && \
      eval "$(micromamba shell hook --shell bash)" && \
      micromamba activate base && \
      micromamba clean --all --yes

RUN eval "$(micromamba shell hook --shell bash)" && \
    micromamba activate base && \
    Rscript -e 'install.packages(c("Matrix", "lme4", "lmerTest"), repos="https://cloud.r-project.org")' && \
    Rscript -e 'devtools::install_github("thiesgehrmann/multidiffabundance@devel", dependencies=FALSE)'
    
LABEL maintainer="Thies Gehrmann"

#ENTRYPOINT ["/usr/bin/mda"] # Do not provide entrypoint for now
