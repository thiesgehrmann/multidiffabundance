FROM mambaorg/micromamba

# Assumes that build environment is in /home/thies/repos/devel/multidiffabundance
COPY ./ /multidiffabundance
COPY ./MDA/mda /use/bin/mda

ENV PATH="/opt/conda/bin:$PATH"

RUN micromamba install \
      -y -n base \
      -c bioconda -c conda-forge \
      r-base \
      r-remotes \
      r-pak \
      r-tidyverse \
      r-digest \
      r-lmerTest \
      r-reshape2  \
      r-vegan \
      r-GUniFrac \
      r-lme4 \
      r-matrix  && \
      eval "$(micromamba shell hook --shell bash)" && \
      micromamba activate base && \
    (echo "remotes::install_github('jsilve24/ALDEx3'); pak::local_install('/multidiffabundance', dependencies=FALSE)" | R --no-save) && \
    micromamba clean --all --yes
    
LABEL maintainer="Thies Gehrmann"

#ENTRYPOINT ["/usr/bin/mda"] # Do not provide entrypoint for now
