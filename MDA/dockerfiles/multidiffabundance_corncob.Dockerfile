FROM mambaorg/micromamba

# Assumes that build environment is in /home/thies/repos/devel/multidiffabundance
COPY ./ /multidiffabundance
COPY ./MDA/mda /use/bin/mda

ENV PATH="/opt/conda/bin:$PATH"

RUN micromamba install \
      -y -n base \
      -c bioconda -c conda-forge \
      cmake \
      r-base \
      r-pak \
      r-tidyverse \
      r-digest \
      r-reshape2  \
      r-vegan \
      r-corncob \
      bioconductor-phyloseq && \
      eval "$(micromamba shell hook --shell bash)" && \
      micromamba activate base && \
      micromamba clean --all --yes

RUN eval "$(micromamba shell hook --shell bash)" && \
    micromamba activate base && \
    Rscript -e 'install.packages(c("Matrix", "lme4", "lmerTest"), repos="https://cloud.r-project.org")' && \
    Rscript -e 'pak::local_install("/multidiffabundance", dependencies=FALSE)'
# Rscript -e 'devtools::install_github("thiesgehrmann/multidiffabundance@devel", dependencies=FALSE)
    
LABEL maintainer="Thies Gehrmann"

#ENTRYPOINT ["/usr/bin/mda"] # Do not provide entrypoint for now
