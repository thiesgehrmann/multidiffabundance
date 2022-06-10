# docker build -t thiesgehrmann/multidifferentialabundance:1 ./
FROM continuumio/miniconda3

ADD MDA /usr/MDA/
RUN chmod +x /usr/MDA/setup_conda_envs.sh && /usr/MDA/setup_conda_envs.sh


LABEL maintainer="Thies Gehrmann"

ENTRYPOINT ["/usr/MDA/mda"]