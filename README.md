# MDA: multidiffabundance
A toolkit for the testing of differential abundance with many different tools, each provided with a similar interface and a compatible output format.
The following packages are currently supported (Only tools that allow for adjustment with other variables are selected): ALDEx2, ANCOMBC, Corncob (optionally), DESeq2, Limma(voom), lm/lmer CLR, Maaslin2 and ZicoSeq.
We also provide functionality to perform alpha and beta diversity tests, and meta-data associations to categorical and continuous variables (typical extracted from the microbiome data)

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
 
Example data is provided in the `mda.example` object provided by the R package, and also in the data files found in `inst/extdata/` of this github repository.

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
`out$res.full` contains the output of all variables per method! `out$res` contains only the chosen variable.

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
 5. `mda.deseq2`: Run DESeq2 online
 6. `mda.limma`: Run Limma only
 7. `mda.lmclr`: Run clr(abundance) ~ model only
 8. `mda.maaslin2`: Run Maaslin2 only
 9. `mda.zicoseq`: Run Zicoseq only

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
    abundance="inst/extdata/mda.example.count_data.tsv"
    meta_data="inst/extdata/mda.example.meta_data.tsv"
    functions="inst/extdata/mda.example.formulas.txt"
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

## In an error message, or a warning, I see a strange variable name that wasn't in my formula. What's up?

mda renames all variables in order to keep track of them throughout all the different tools. The order of the original formula is preserved, but all terms in your formula are expanded, so the number of variables may change (e.g. with interaction effects or with categorical variables). Each term in your formula is given a new name. You can look up the name in, for example, the first formula in your object: `D$formula[[1]]$map`.

## How can I limit the number of CPUs used?
Some tools are parallelized and use all the CPUs available - even not making it possible to change this with a setting. However, `taskset` comes to the rescue!
Run the mda command as so:

```shell
    taskset -c 0-8 ./mda blah blah blah
```

## Which tools are able to adjust for covariates, and which are able to accept random effects?

While all models can accept adjustment terms ala a linear model, only some tools can accept random intercept effects or random slopes. If you attempt to evaluate a formula with random effects using a tool that cannot handle it, MDA returns a dummy output. For these tools, if you want to adjust for a repeated measures, you may need to encode it as a fixed effect in order to make use of all the tools (or better: don't use the tools).

| Tool          | Covariates             | Random intercepts      | Random slopes          | Effects reported       |
| ------------- |:----------------------:|:----------------------:|:----------------------:|:----------------------:|
| ALDEx2        | ✓                      |                        |                        | All                    |
| ANCOM-BC      | ✓                      | ✓ (many)               | ✓ (many)               | All                    |
| Corncob       | ✓                      |                        |                        | All                    |
| DESeq2        | ✓                      |                        |                        | All                    |
| limma         | ✓                      | ✓ (one)                |                        | All                    |
| lmCLR         | ✓                      | ✓ (many)               | ✓ (many)               | All                    |
| Maaslin2      | ✓                      | ✓ (many)               |                        | All                    |
| Zicoseq       | ✓                      |                        |                        | First                  |
| alpha         | ✓                      | ✓ (many)               | ✓ (many)               | All                    |
| beta          | ✓                      | ✓ (one) \*see note     |                        | First \*\*see note     |
| continuous    | ✓                      | ✓ (many)               | ✓ (many)               | All                    |
| group         | ✓                      | ✓ (many)               | ✓ (many)               | All                    |

\* In adonis2, this is not a real random intercept. It merely performs permutations within a certain grouping as defined by a categorical variable. This can be interpreted as considering this grouping as a random effect.

\*\* This is due to a hack implemented to speed up the test. If, when running mda.beta, you set `beta.hack=FALSE`, then it will take longer but return all tested effects.


In the future, I am hoping also to include `dream` in this list, which is essentially limma, but allows multiple random intercepts.

## How should I format my variables of interest

mda does exensive processing of formulas.
Essentially, you are able to provide any formula, no matter how excessively complex. e.g.

```R
~ log2(abs(I(a*b+c))) * z + (d+e+f)^2 + (1|q) + (exp(m)|g)
```

However, different tools handle these formulas differently, and some do not provide results for all formulas.
mda returns a list with two dataframes:
 * `res` : A table in which **ONLY THE RESULTS FOR THE FIRST VARIABLE OF THE FORMULA ARE PROVIDED**
 * `res.full`: A table in which all results for all variables (if provided by the tool) are provided
 
If, for example, you use a function like: `~ a*b`, this is expanded to `~ a + b + a:b`, in this order.
(Internally, these are then renamed to ` mda_00001 + mda_00002 + mda_00003`)
If your effect of interest is `a:b`, then this effect will not be reported in `res`
Therefore, the `res` results may not be what you are looking for. Inspect the `res.full` dataframe.
Alternatively, you can specify the order manually, by specifying this in the formula: `~ a:b + a + b`.

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

## Can you include my favourite method `xxx`?
Maybe. Send me a message.

## How do I deal with all the results?

Good question! Our current approach is to count the number of tools that report a significant result. This gives an indication of the consensus of the tools. For example, when 5/6 tools agree on a significant effect, I am more confident in the result than if only 2/6 tools agree. Our threshold for reporting is at least 3/6 tools.

## Can you include method `xxx`?
Maybe. Send me a message.
