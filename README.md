# MDA: multidiffabundance
A toolkit for the testing of differential abundance with many different tools, each provided with a similar interface and a compatible output format.

# Installation

## Installation of dependencies

 MDA has many dependencies (as it is a collection of so many tools).
 They (WILL BE) automatically installed when installing the package via github as below.
 However, you can also install the dependencies with conda:
 
```bash
  wget https://github.com/thiesgehrmann/multidiffabundance/blob/main/MDA/conda_deps.yml ./mda_conda_deps.yml
  conda env create -n mda --file ./mda_conda_deps.yml
```

## Installation of R package

 MDA has the following 

 Install the R package with devtools
 
 ```R
    # install.packages("devtools")
    devtools::install_github("thiesgehrmann/multidiffabundance", dependencies=TRUE)
 ```
 
## Installation of the command line tool

 To use the command line tool, you should have the R package installed.
 Then run the following commands:
 
```bash
  wget https://raw.githubusercontent.com/thiesgehrmann/multidiffabundance/main/MDA/mda ./
  chmod +x ./mda
  sudo mv ./mda /usr/bin # not necessary
```

## Installation of the Docker/Singularity image

 We provide a wrapper for the docker image, in which all necessary dependencies are installed
 For this, you do not need to install the package

# Running the tool

The input to MDA is the following:
 1. A table of samples x taxa (rows x columns)
 2. A table of samples x metadata (rows x columns)
 3. A formula, or a file with a list of formulas
 4. A path to an output folder


## Running in R

(The interface to these functions in R will be improved in future versions)

```R
    args <- c("data/moving-pics-abundance-cleaned.tsv",
              "data/moving-pics-meta-cleaned.tsv",
              "data/moving-pics-functions.txt",
              "output_folder")
    D <- mda.load(args)
    out <- mda.all(D)
    out$res # Relevant output data here
```

## Running via the command line

```bash
    # Assumes R package is installed
    abundance="data/moving-pics-abundance-cleaned.tsv"
    meta_data="data/moving-pics-meta-cleaned.tsv"
    functions="data/moving-pics-functions.txt"
    outdir="output_folder"
    
    wget https://raw.githubusercontent.com/thiesgehrmann/multidiffabundance/main/MDA/mda ./
    chmod +x ./mda
    ./mda "$abundance" "$meta_data" "$functions" "$outdir" # Produces output in $outdir/results.tsv
```

## Running via docker/singularity image

```bash
    # Assumes Docker or singularity is installed on your compute node
    abundance="data/moving-pics-abundance-cleaned.tsv"
    meta_data="data/moving-pics-meta-cleaned.tsv"
    functions="data/moving-pics-functions.txt"
    outdir="output_folder"
    docker_mda "$abundance" "$meta_data" "$functions" "$outdir" # Produces output in $outdir/results.tsv
    
```
