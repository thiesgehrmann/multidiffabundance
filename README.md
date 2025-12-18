# MDA: multidiffabundance
A toolkit for the testing of differential abundance with many different tools, each provided with a similar interface and a compatible output format.
The following packages are currently supported (Only tools that allow for adjustment with other variables are selected): ALDEx2, ANCOMBC, Corncob, DESeq2, Limma(voom), lm/lmer CLR, and Maaslin2.

# Quickstart

```R
# install.packages("devtools")
devtools::install_github("thiesgehrmann/multidiffabundance", dependencies=TRUE)

library(multidiffabundance)
data("mda.example", package="multidiffabundance")
D <- mda.create(mda.example$count_data, # A dataframe of sample taxa counts (samples are rows)
                mda.example$meta_data,  # A dataframe of sample meta data (samples are rows, rowname should match count data dataframe)
                mda.example$formulas,   # A list of formulas (or a single formula)
                usecache=TRUE,          # Cache computationally expensive steps (default behaviour)
                recache=FALSE)          # Do not overwrite existing cache files (default behaviour)
out <- mda.all(D) # Runs all methods (you can specify which you want to run)
out$res           # Relevant output data here
```


# Installation

## Installation of dependencies

 MDA has many dependencies (as it is a collection of so many tools).
 They (SHOULD BE) automatically installed when installing the package via github as in the section below.
 However, you can also install the dependencies with conda:
 
```bash
wget https://raw.githubusercontent.com/thiesgehrmann/multidiffabundance/refs/heads/main/MDA/mda_conda_deps.yml
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
    mda 
```

## Installation of the Docker/Singularity image

 We provide a wrapper for the docker image, in which all necessary dependencies are installed
 For this, you do not need to install anything (other than docker or singularity).
 
```shell
    wget https://raw.githubusercontent.com/thiesgehrmann/multidiffabundance/main/MDA/container_mda.sh
    chmod +x ./container_mda.sh
    sudo mv ./container_mda.sh /usr/bin # not necessary
```

# Running the tool

The input to MDA is the following:
 1. A table of samples x taxa (rows x columns, first column should be sample ID)
 2. A table of samples x metadata (rows x columns, first column should be sample ID)
    Column names should not be 
 3. A (list of) formula(s), or a file with a line-separated list of formulas
 4. A path to an output folder (optional except in command line version)
 
Example data is provided in the `mda.example` object provided by the R package, and also in the data files found in `data/` of this github repository.

## Running in R

```R
    library(multidiffabundance)
    data("mda.example", package="multidiffabundance")
    D <- mda.create(mda.example$count_data,
                    mda.example$meta_data,
                    mda.example$formulas)
    out <- mda.all(D)
    out$res # Relevant output data here
```

**The effect for the FIRST variable in the formula will be reported**, so make sure to format your formulas correctly.
`out$res.full` contains the output of all variables per method, but it is not cleaned! `out$res` contains only the chosen variable.

Functions to create the MDA objects are:
 1. `mda.create(counts, meta, formulas, ...)`: Takes loaded objects
 2. `mda.from_cmdargs(c(counts, meta, formulas))` : Takes command line arguments
 3. `mda.from_files(counts_file, meta_file, formulas, ...)` : Takes files
 3. `mda.from_tidyamplicons(ta, formulas, ...)` : Takes a tidyamplicons object
 4. `mda.from_phyloseq(ps, formulas, ...)` : Takes a phyloseq object (unimplemented)
 
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
    abundance="data/mda.example.count_data.tsv" # The following files exist in the github repository under data/
    meta_data="data/mda.example.meta_data.tsv"
    functions="data/mda.example.formulas.txt"
    outdir="output_folder"
    
    mda "$abundance" \
        "$meta_data" \
        "$functions" \
        "$outdir" | tee "$outdir/stdout" # Produces output in "$outdir/results.tsv"
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

# Frequently asked questions?

## How can I limit the number of CPUs used?
Some tools are parallelized and use all the CPUs available - even not making it possible to change this with a setting. However, `taskset` comes to the rescue!
Run the mda command as so:

```shell
    taskset -c 0-8 ./mda blah blah blah
```

## Which tools are able to adjust for covariates, and which are able to accept random effects?

While all models can accept adjustment terms ala a linear model, only three tools can accept random intercept effects, and only one can accept random slopes. If you attempt to evaluate a formula with random effects using a tool that cannot handle it, MDA issues a warning but DOES NOT FAIL. It evaluates the formula without the random terms. limma can only handle one random intercept effect. If you want to adjust for a repeated measures, you may need to encode it as a fixed effect in order to make use of all the tools.

| Tool          | Covariates             | Random intercepts      | Random slopes          |
| ------------- |:----------------------:|:----------------------:|:----------------------:|
| ALDEx2        | ✓                      |                        |                        |
| ANCOM-BC      | ✓                      |                        |                        |
| Corncob       | ✓                      |                        |                        |
| DESeq2        | ✓                      |                        |                        |
| limma         | ✓                      | ✓ (one)                |                        |
| lmCLR         | ✓                      | ✓ (many)               | ✓ (many)               |
| Maaslin2      | ✓                      | ✓ (many)               |                        |

In the future, I am hoping also to include `dream` in this list, which allows multiple random intercept

## How should I format my variables of interest

In this implementation (which is maybe not the best, and perhaps it will change), the first variable of each provided formula is the variable of interest, and only the effect of this variable will be reported.
In the future there may be ways to retrieve multiple effects.
For now, to retrieve the effects of multiple variables from the same formula, you must provide multiple formulas in which only the order of the variables has been changed. You can use the `mda.permute_formula` function for a crude implementation of this approach.
```R
mda.permute_formula(~a+b+c) # results in ~a+b+c, ~b+a+c, ~c+a+b
```

This should not result in increased computational overhead as the fitted models are cached with an formula-term-order-invariant hash (i.e. `~a+b` will be stored in the same cache as `~b+a`)

```R
    form <- ~ ReportedAntibioticUsage + DaysSinceExperimentStart + (1|Subject)
    D <- mda.create(
             mda.example$count_data,
             mda.example$meta_data,
             mda.permute_formula(form))

    mda.maaslin2(D)
```

A correlary of the previous statement is that if you are evaluating categorical variables, only one of the dummy encoded results will be reported.
It may therefore be wise to dummycode your variables in advance in order to maintain control over the contrasts. Try out the `fastDummies` package for this.

## Why is it a problem if a variable has a name that is a prefix of another?

Currently, to extract the variable of interest, we select it from the output of the method. In the case of categorical variables, this can be tricky. Say we have a categorical variable called `var1` with values `TRUE` and `FALSE`. It is usually listed in the output as `var1TRUE`. Therefore, we select the variable that starts with `var1`. This can be a problem if you have multiple variables and one of them has your variable of interest as a prefix. For example, your variable of interest is `var1`, but you also have `var12`. In this case, we may accidentally select `var12` instead of `var1`. In this case, please rename your variable so that it is not the prefix of another. This issue may be solved later.

## How can I modify the caching behaviour?

By default, the computationally expensive model runs are cached in a way that is explicitly linked to the method used, the count data, the meta data, and the formula currently being evaluated. You can turn off model caching with the `--nocache` option on the commandline, or by setting `usecache=FALSE` in the `mda.create` functions.

You can also overwrite existing cache files with `--recache` option on the commandline, or by setting `recache=TRUE` in the `mda.create` functions.

## I'd like to modify the default parameters to some of the tools used - how can I do that.
Currently this is not implemented. At some point I want to store the default parameters in an R object that can be modified in the mda data object, but not yet. For now you will need to modify the source code. This is easily done (though I agree it is a pain):

```bash
git clone https://github.com/thiesgehrmann/multidiffabundance.git
cd multidiffabundance
# *make changes to R/mda.method.R*
echo "devtools:install_local(force=TRUE, dependencies=FALSE)" | R --no-save
```

## 

## Can you include method `xxx`?
Maybe. Send me a message.
