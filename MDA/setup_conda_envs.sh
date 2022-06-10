#!/usr/bin/env bash
# Script to set up the conda environments for the MDA tool

conda env create --file /usr/MDA/MDA/aldex2/conda_deps.yml # For ALDEX2
conda env create --file /usr/MDA/MDA/ancombc/conda_deps.yml # For ANCOMBC
conda env create --file /usr/MDA/MDA/corncob/conda_deps.yml # For CORNCOB
conda env create --file /usr/MDA/MDA/deseq2/conda_deps.yml # For DESEQ
conda env create --file /usr/MDA/MDA/limma/conda_deps.yml # For LIMMA
conda env create --file /usr/MDA/MDA/limma/conda_deps.yml # For LMCLR, is same as for limma
conda env create --file /usr/MDA/MDA/maaslin2/conda_deps.yml # For MAASLIN2
conda clean