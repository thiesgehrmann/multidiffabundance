# MDA: multidiffabundance
A toolkit for the testing of differential abundance with many different tools, each provided with a similar interface and a compatible output format.
The following packages are currently supported (Only tools that allow for adjustment with other variables are selected): ALDEx2, ANCOMBC, Corncob, DESeq2, Limma(voom), lm CLR, and Maaslin2.

# Quickstart

```R
# install.packages("devtools")
devtools::install_github("thiesgehrmann/multidiffabundance", dependencies=TRUE)

library(multidiffabundance)
data("mda.example", package="multidiffabundance")
D <- mda.create(mda.example$count_data,
                mda.example$meta_data,
                mda.example$formulas)
out <- mda.all(D) # Runs all methods
out$res # Relevant output data here
```


# Installation

## Installation of dependencies

 MDA has many dependencies (as it is a collection of so many tools).
 They (SHOULD BE) automatically installed when installing the package via github as in the section below.
 However, you can also install the dependencies with conda:
 
```bash
wget https://github.com/thiesgehrmann/multidiffabundance/blob/main/MDA/mda_conda_deps.yml ./mda_conda_deps.yml
conda env create -n mda --file ./mda_conda_deps.yml
```

However, it is also possible to install MDA without all the dependencies, and only use the packages that you have installed.

## Installation of R package

 Install the R package with devtools:
 
 ```R
# install.packages("devtools")
devtools::install_github("thiesgehrmann/multidiffabundance", dependencies=TRUE)
 ```
 
## Installation of the command line tool

 To use the command line tool, you should have the R package installed.
 Then run the following commands:
 
```shell
wget https://raw.githubusercontent.com/thiesgehrmann/multidiffabundance/main/MDA/mda ./
chmod +x ./mda
sudo mv ./mda /usr/bin # not necessary
```

## Installation of the Docker/Singularity image

 We provide a wrapper for the docker image, in which all necessary dependencies are installed
 For this, you do not need to install anything (other than docker or singularity)
 
```shell
wget https://raw.githubusercontent.com/thiesgehrmann/multidiffabundance/main/MDA/docker_mda.sh
chmod +x ./docker_mda.sh
sudo mv ./docker_mda.sh /usr/bin # not necessary
```

# Running the tool

The input to MDA is the following:
 1. A table of samples x taxa (rows x columns, first column should be sample ID)
 2. A table of samples x metadata (rows x columns, first column should be sample ID)
 3. A (list of) formula(s), or a file with a line-separated list of formulas
 4. A path to an output folder (optional except in command line version)


## Running in R

(The interface to these functions in R will be improved in future versions)

```R
library(multidiffabundance)
data("mda.example", package="multidiffabundance")
D <- mda.create(mda.example$count_data,
                mda.example$meta_data,
                mda.example$formulas)
out <- mda.all(D)
out$res # Relevant output data here
```

**The effect for the FIRST variable in the formula will be reported**, so make sure to format your formulas. `out$res.all` contains the output of all variables per method, but it is not cleaned!

Functions to create the MDA objects are:
 1. `mda.create`: Takes loaded objects
 2. `mda.from_cmdargs` : Takes command line arguments
 3. `mda.from_files` : Takes files
 3. `mda.from_tidyamplicons` : Takes a tidyamplicons object
 4. `mda.from_phyloseq` : Takes a phyloseq object (unimplemented)
 
Functions to run the differential abundance tests are:
 1. `mda.all`: Run all tools, or a combination of tools
 2. `mda.aldex2`: Run ALDEx2 only
 3. `mda.ancom`: Run ANCOMBC only
 4. `mda.corncob`: Run Corncob only
 5. `mra.deseq2`: Run DESeq2 online
 6. `mra.limma`: Run Limma only
 7. `mra.lmclr`: Run clr(abundance) ~ model only
 8. `mra.maaslin2`: Run Maaslin2 only

## Running via the command line

```shell
# Assumes R package is installed
abundance="data/mda.example.count_data.tsv"
meta_data="data/mda.example.meta_data.tsv"
functions="data/mda.example.formulas.txt"
outdir="output_folder"

mda "$abundance" \
      "$meta_data" \
      "$functions" \
      "$outdir" # Produces output in $outdir/results.tsv
```

## Running via docker/singularity image

A docker image is [provided on dockerhub](https://hub.docker.com/repository/docker/thiesgehrmann/multidiffabundance), as defined by the [Dockerfile](https://raw.githubusercontent.com/thiesgehrmann/multidiffabundance/main/MDA/Dockerfile).

A wrapper script `container_mda.sh` allows you to run the docker image with the same interface as the standalone script above. The script uses singularity by default, but you can also specify to use docker with the --docker argument.

```shell
# Assumes Docker or singularity is installed
abundance="data/mda.example.count_data.tsv"
meta_data="data/mda.example.meta_data.tsv"
functions="data/mda.example.formulas.txt"
outdir="output_folder"

# Runs with singularity
container_mda.sh \
    "$abundance" \
    "$meta_data" \
    "$functions" \
    "$outdir" # Produces output in $outdir/results.tsv

# Runs with docker
container_mda.sh --docker \
    "$abundance" \
    "$meta_data" \
    "$functions" \
    "$outdir" # Produces output in $outdir/results.tsv

```
