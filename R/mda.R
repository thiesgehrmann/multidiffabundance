###############################################################################
# MDA Common R functions
###############################################################################
                                                      
###############################################################################
# Input loading functions

mda.from_cmdargs <- function(args, ...){
    require(tools)
    require("tidyverse")
    
    if (length(args) != 4){
        message("[MDA] mda.from_cmdargs ERROR: Too few arguments. There should be 4!")
        quit(1)
    }
    
    abundance <- args[1]
    meta <- args[2]
    formula.data <- args[3]
    outprefix <- args[4]

    mda.from_files(abundance, meta, formula.data, outprefix)
}
                                                      
mda.from_files <- function(abundance, meta, formula.data, outprefix=tempdir(), ...){
    # Prepare formula
    require(tidyverse)

    raw.formula.data <- mda.load_formula_input(formula.data)
    if (!mda.verify_formula_input(raw.formula.data)){
        message("[MDA] mda.from_files ERROR: Errors in processing formula input.")
        quit(1)
    }

    ###############################################################################
    # Load the data

    count_data <- as.data.frame(read_tsv(abundance))
    row.names(count_data) <- count_data[,1]
    count_data <- count_data[,-1]

    meta_data <- as.data.frame(read_tsv(meta))
    row.names(meta_data) <- meta_data[,1]
    meta_data <- meta_data[,-1]

    common_samples <- intersect(rownames(count_data), rownames(meta_data))
    meta_data <- meta_data[common_samples, ]
                                                        
    ###############################################################################
    
    dat <- mda.create(count_data, meta_data, raw.formula.data, outprefix, ...)
    
    return(dat)
}

mda.from_tidyamplicons <- function(ta, formulas, output_dir=tempdir(), ...){
    require(dpylr)
    require(tidyr)
    meta_data <- as.data.frame(ta$samples)
    rownames(meta_data) <- meta_data$sample_id
    
    count_data <- as.data.frame(pivot_wider(ta$abundances, id_cols="sample_id", names_from="taxon_id", values_from="abundance", values_fill=0))
    rownames(count_data) <- count_data$sample_id
    count_data <- count_data[,-1]
    
    dat <- mda.create(count_data, meta_data, formulas, output_dir, ...)
    dat$ta <- ta
    
    return(dat)
}

mda.from_phyloseq <- function(phys, formulas, output_dir=tempdir(), ...){
    message("[MDA] mda.from_phyloseq ERROR: UNIMPLEMENTED")
}
                                                      
mda.create <- function(count_data, meta_data, formulas, output_dir=tempdir(), usecache=TRUE, recache=FALSE){
    require(digest)
    if (! all(rownames(count_data) == rownames(meta_data))) {
        message("[MDA] mda.create ERROR: Rownames of meta data and count data do not match!")
        return(NULL)
    }
    
    FD <- mda.process_formula_input(unlist(lapply(c(formulas), function(x){mda.deparse(x)})))
    nonrare <- mda.nonrare_taxa(count_data, 0.1) 

    numeric_meta <- colnames(meta_data)[unlist(lapply(colnames(meta_data), function(x)is.numeric(meta_data[,x])))]
    meta_data[,numeric_meta] <- scale(meta_data[,numeric_meta])
    
    dat <- c()
    dat$count_data  <- count_data
    dat$nonrare     <- nonrare
    dat$meta_data   <- meta_data
    dat$formula     <- FD
    dat$outprefix   <- output_dir
    dat$usecache    <- usecache
    dat$recache     <- recache
    checksums <- unlist(lapply(list(count_data, meta_data), digest, algo="md5"))
    dat$cacheprefix <- paste0(c(output_dir, "mda.cache", paste0(checksums, collapse=".")), collapse="/")
    if (usecache){
        mda.mkdirp(dirname(dat$cacheprefix))
    }

    return(dat)
}

###############################################################################
# Filesystem functionality

mda.mkdirp <- function(dir){
    if (!dir.exists(dir)){
        mda.mkdirp(dirname(dir))
        dir.create(dir)
    }
}

###############################################################################
# Formula processing functions
                                                      
mda.load_formula_input <- function(formula_input){
  raw <- if (inherits(formula_input, "formula")) {
      mda.deparse(formula_input)
  } else if (file.exists(formula_input)){
    unlist(strsplit(readLines(formula_input), "\n"))
  } else {
    formula_input
  }
  raw
}

mda.verify_formula_input <- function(raw_formula){
  
  valids <- lapply(raw_formula, function(rs){
    tryCatch({as.formula(rs); TRUE},
             warning=function(w){FALSE},
             error=function(r){FALSE})
  })
  all(unlist(valids))
}

mda.process_formula_input <- function(raw_formula){
  form <- lapply(raw_formula, as.formula)
  
  fterms <-                  lapply(form, function(x){labels(terms(x))})
  fterms_fixed  <-           lapply(fterms, function(x){ x[!grepl("\\|", x)]})
  fterms_random_intercept <- lapply(fterms, function(x){ x[grepl("^[1][ ]*\\|", x)]})
  fterms_random_slope <-     lapply(fterms, function(x){ x[grepl("^[01]?[ ]*[+]?[^0-9].+\\|", x)]})

  FD <- c()
  FD$raw <- raw_formula
  FD$formula <- form
  FD$main_var <- unlist(lapply(fterms_fixed, function(x){unlist(x[1])}))
  FD$adj_vars <- lapply(fterms_fixed, function(x){paste0(x[-1], collapse="+")})

  FD$rand_intercept <- lapply(fterms_random_intercept, function(ftr){
    unlist(lapply(ftr, function(x){paste0(c('(', x, ')'), collapse="")}))
  })
    
  FD$rand_slope <- lapply(fterms_random_slope, function(ftr){
    unlist(lapply(ftr, function(x){paste0(c('(', x, ')'), collapse="")}))
  })
    
  FD$norand <- lapply(fterms_fixed, function(ftf){
    as.formula(paste0(c("~",paste0(ftf, collapse="+")),collapse=""))
  })
  
  #if (nchar(paste0(unlist(fterms_random), collapse="")) > 0){
  #  message("[MDA] mda.process_formula_input WARNING: No mixed effect terms are allowed in this implementation. Random effects will be ignored. Please format your random effects as fixed effects. Sorry.")
  #}
  
  FD
}
                                                      
mda.permute_formula <- function(form){
    unlist(lapply(c(form), function(f){
        L <- labels(terms(as.formula(f)))
        Lf <- L[!grepl("\\|",L)]
        Lr <- L[grepl("\\|",L)]
        Lrf <- unlist(lapply(Lr, function(x) paste0(c("(",x,")"), collapse="")))
        unlist(lapply(Lf, function(l){
            other <- setdiff(Lf, l)
            as.formula(paste0(c(paste0(c("~",l,""), collapse=""), other, Lrf), collapse="+"))
        }))
    }))
}

###############################################################################
# Abundance processing functions


mda.nonrare_taxa <- function(table , cutoff_pct) {
    cutoff  <- ceiling(cutoff_pct * nrow(table))
    nonzero <- colnames(table)[colSums(table > 0) >= cutoff]
    return(nonzero)
}
                                                      
mda.relative_abundance <- function(count_data){
    depth <- rowSums(count_data)
    radata <- apply(count_data, 2, function(x){x/depth})
    return(radata)
}


mda.pseudocount <- function(count_data){
    pcount <- (rowSums(count_data) / max(rowSums(count_data)))
    pdata <- apply(count_data,2,function(x){x+pcount})
    return(pdata)
}
                                                      
mda.clr <- function(df){
    denom <- exp(rowMeans(log(df)))
    df / denom
}

###############################################################################
# Statistical counts functions
                                                      
mda.meta.n <- function(D, var){
    sum(unlist(lapply(D$meta_data[,var],function(x){!(is.na(x))})))
}


mda.meta.freq <- function(D, var){
    if(is.numeric(D$meta_data[,var])){
        ""
    } else {
        paste0(mapply(function(x,y){paste0(c(x,y), collapse=":")}, names(table(D$meta_data[,var])), as.character(table(D$meta_data[,var]))), collapse=", ")
    }
}
                                                      
###############################################################################
# output cache functions
                                                      
mda.cache_filename <- function(outprefix, method, form, suffix="tsv", collapse="."){
    require(digest)
    mainvar <- labels(terms(as.formula(form)))[1]
    L <- labels(terms(as.formula(form)))
    form.fmt <- paste0(sort(L), collapse="+")
    
    form.digest <- digest(form.fmt, algo="md5")
    form.digest <- gsub("[~():._! |]", "", form.digest)
    
    filename <- paste0(c(paste0(c(outprefix, method, form.digest), collapse=collapse), suffix), collapse=".")
    filename
}

mda.cache_save <- function(dat, outprefix, method, form, suffix="rds", ...){
  filename <- mda.get_cache_filename(outprefix, method, form, ...)
  saveRDS(dat, filename)
}

mda.cache_load <- function(outprefix, method, form, suffix="rds"){
  filename <- mda.get_cache_filename(outprefix, method, form, ...)
  dat <- readRDS(filename)
}
                                                      
mda.cache_load_or_run_save <- function(mda.D, method, form, expr) {
    D <- mda.D
    outputprefix <- D$cacheprefix
    cache.file <- mda.cache_filename(outputprefix, method, form, suffix="rds")
    mainvar <- labels(terms(as.formula(form)))[1]
    data <- if (file.exists(cache.file) & (D$usecache) & !(D$recache)){
        message(paste0(c("[MDA] CacheLoad: ", method, ", ", mainvar, " (", basename(cache.file), ")"), collapse=""))
        readRDS(cache.file)
    } else{
        data <- expr
        if (D$usecache){
            message(paste0(c("[MDA] CacheStore: ", method, ", ", mainvar, " (", basename(cache.file), ")"), collapse=""))
            saveRDS(data, cache.file)
        }
        data
    }
    data
}
                                                      
mda.deparse <- function(form){
    if (inherits(form, what="formula")){
        # Taken from as.character.formula in formula.tools
        form <- paste( deparse(form), collapse=" " )
        form <- gsub( "\\s+", " ", form, perl=FALSE ) # remove multiple spaces
        return(form)
    } else { # Assume it is a string already
        return(form)
    }
}
                                                      
###############################################################################