# docker build -t thiesgehrmann/multidiffabundance:1 ./
# docker push thiesgehrmann/multidiffabundance:1
FROM continuumio/miniconda3

COPY mda /usr/bin/mda
COPY mda_conda_deps.yml /usr/mda_conda_deps.yml
RUN chmod +x /usr/bin/mda && \
    conda update -n base -c defaults conda && \
    conda env update -n base --file /usr/mda_conda_deps.yml && \
    rm /usr/mda_conda_deps.yml && \
    (echo "library(devtools); devtools::install_github('thiesgehrmann/multidiffabundance', dependencies=FALSE)" | R --no-save) && \
    conda clean -a

LABEL maintainer="Thies Gehrmann"

#ENTRYPOINT ["/usr/bin/mda"] # Do not provide entrypoint for now