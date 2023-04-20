###############################################################################
# MDA Common R functions
###############################################################################
                                                      
###############################################################################
# Input loading functions

mda.from_cmdargs <- function(args, ...){
    suppressPackageStartupMessages(require(tools))
    suppressPackageStartupMessages(require("tidyverse"))
    
    if (length(args) != 4){
        message("[MDA] mda.from_cmdargs ERROR: Too few arguments. There should be 4!")
        quit(1)
    }
    
    abundance <- args[1]
    meta <- args[2]
    formula.data <- args[3]
    outprefix <- args[4]

    mda.from_files(abundance, meta, formula.data, outprefix, ...)
}
                                                      
mda.from_files <- function(abundance, meta, formula.data, outprefix=tempdir(), ...){
    # Prepare formula
    suppressPackageStartupMessages(require(tidyverse))

    raw.formula.data <- mda.load_formula_input(formula.data)
    if (!mda.verify_formula_input(raw.formula.data)){
        message("[MDA] mda.from_files ERROR: Errors in processing formula input.")
        quit(1)
    }

    ###############################################################################
    # Load the data

    count_data <- as.data.frame(read_tsv(abundance, col_types = cols()))
    row.names(count_data) <- count_data[,1]
    count_data <- count_data[,-1]

    meta_data <- as.data.frame(read_tsv(meta, col_types = cols()))
    row.names(meta_data) <- meta_data[,1]
    meta_data <- meta_data[,-1]

    common_samples <- intersect(rownames(count_data), rownames(meta_data))
    meta_data <- meta_data[common_samples, ]
                                                        
    ###############################################################################
    
    dat <- mda.create(count_data, meta_data, raw.formula.data, outprefix, ...)
    
    return(dat)
}

mda.from_tidyamplicons <- function(ta, formulas, output_dir=tempdir(), ...){
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tidyr))
    meta_data <- as.data.frame(ta$samples)
    rownames(meta_data) <- meta_data$sample_id
    
    count_data <- as.data.frame(pivot_wider(ta$abundances, id_cols="sample_id", names_from="taxon_id", values_from="abundance", values_fill=0))
    rownames(count_data) <- count_data$sample_id
    count_data <- count_data[,-1]

    union.samples = sort(union(rownames(count_data), rownames(meta_data)))
    intersect.samples = sort(intersect(rownames(count_data), rownames(meta_data)))
    
    if ( length(union.samples) != length(intersect.samples) ){
        message("[MDA] mda.from_tidyamplicons WARNING: tidyamplicons object has differing samples present in abundance and meta data!")
    }
    
    dat <- mda.create(count_data[intersect.samples,], meta_data[intersect.samples,], formulas, output_dir, ...)
    dat$ta <- ta
    
    return(dat)
}

mda.from_phyloseq <- function(phys, formulas, output_dir=tempdir(), ...){
    message("[MDA] mda.from_phyloseq ERROR: UNIMPLEMENTED")
}
                                                      
mda.create <- function(count_data, meta_data, formulas, output_dir=tempdir(), usecache=TRUE, recache=FALSE, nonrare.pct=0.1){
    suppressPackageStartupMessages(require(digest))
    
    if (! all(rownames(count_data) == rownames(meta_data))) {
        message("[MDA] mda.create ERROR: Rownames of meta data and count data do not match!")
        return(NULL)
    }
    
    FD <- mda.process_formula_input(unlist(lapply(list(formulas), function(x){mda.deparse(x)})))
    nonrare <- mda.nonrare_taxa(count_data, nonrare.pct) 

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
    checksum <- digest(c(count_data, meta_data), algo="md5")
    dat$cacheprefix <- paste0(c(output_dir, "mda.cache", paste0(c(checksum, as.character(nonrare.pct)), collapse="-")), collapse="/")
    if (usecache){
        mda.mkdirp(dirname(dat$cacheprefix))
    }

    return(dat)
}
                                                      
###############################################################################
# Export functions

mda.as_phyloseq <- function(D){
    OTU <- phyloseq::otu_table(t(D$count_data), taxa_are_rows = T)
    sampledata <- phyloseq::sample_data(D$meta_data, errorIfNULL = F)
    phylo <- phyloseq::merge_phyloseq(OTU, sampledata)
    phylo
}
                                                      
mda.as_tidyamplicons <- function(D){
    T
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
                                                      
mda.formula_variables <- function(fn, depth=0, maxdepth=10){
    # Extracts all the variables from a formula, from arithmetic and random effects
    # fn <- ~ log(a) + b*c + (1|dn) + (e|f*g) + I(h+i) + log(j+I(k+l))
    # mda.formula_variables(fn)
    # >> c('a','b','c','dn','e','f','g','h','i','j','k','l')

    if (depth > maxdepth){
        return(NULL)
    }
    vars <- as.character(attributes((terms(as.formula(fn))))$variables)[-1]
    #message(depth, vars)
    subvars <- if (length(vars) == 1){
        vs_random <- unlist(strsplit(vars, '[|]'))
        if (length(vs_random) > 1){
            if (vs_random[1] == "1"){
                subfunc <- as.formula(paste0(c('~ ', vs_random[1]), collapse=""))
                return(mda.formula_variables(subfunc, depth+1,maxdepth))
            } else {
                subfunc <- as.formula(paste0(c('~ ', vs_random[1]," + ",vs_random[2]), collapse=""))
                #message(depth, "paren", subfunc)
                return(mda.formula_variables(subfunc, depth+1,maxdepth))
            }
        }
        
        vs_func <- unlist(strsplit(vars, '[(]'))
        if (length(vs_func) > 1){
            stripfunc <- paste0(vs_func[-1], collapse='(')
            stripfunc <- substr(stripfunc, 1, nchar(stripfunc)-1)
            subfunc <- as.formula(paste0(c('~ ',stripfunc), collapse=""))
            #message(depth, 'paren', subfunc)
            return(mda.formula_variables(subfunc,depth+1,maxdepth))
        }
        vs_random <- unlist(strsplit(vars, '[|]'))

        return(vars)
    } else {
        unlist(lapply(vars, function(x){
            subfunc <- as.formula(paste0(c('~ ',x), collapse=""))
            mda.formula_variables(subfunc,depth+1,maxdepth)
        }))
    }
    subvars
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
  FD$variables <- lapply(FD$raw, mda.formula_variables)

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
                                                      
mda.cache_filename <- function(mda.D, method, form, suffix="tsv", collapse=".", order_invariant=TRUE){
    D <- mda.D
    suppressPackageStartupMessages(require(digest))
    mainvar <- labels(terms(as.formula(form)))[1]
    L <- labels(terms(as.formula(form)))
    L <- if (order_invariant){sort(L)} else{L}
    form.fmt <- paste0(L, collapse="+")
    
    form.digest <- gsub("[~():._! |]", "", form.fmt)
    form.digest <- digest(form.digest, algo="md5")
    data.digest <- digest(c(D$count_data, D$meta_data), algo="md5")
    
    filename <- paste0(c(paste0(c(D$cacheprefix, method, form.digest, data.digest), collapse='_'), suffix), collapse=".")
    filename
}
                                                      
mda.cache_load_or_run_save <- function(mda.D, method, form, expr, order_invariant=TRUE) {
    D <- mda.D
    cache.file <- mda.cache_filename(D, method, form, suffix="rds", order_invariant=order_invariant)
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
# Summary function

mda.summary <- function(res, id_cols="taxa", names_from="variable", method_from="method", pvalue_from="pvalue", qvalue_from="qvalue", effectsize_from="effectsize",
                        values_fn=list, pvalue_threshold=0.05, qvalue_threshold=0.05){

    res <- arrange(res, "method", "taxa", "variable")
    
    list(
        nsig =       pivot_wider(res, id_cols=id_cols, names_from=names_from, values_from=all_of(qvalue_from),     values_fn=function(v){as.integer(sum(v < qvalue_threshold))}),
        pvalue =     pivot_wider(res, id_cols=id_cols, names_from=names_from, values_from=all_of(pvalue_from),     values_fn=list),
        qvalue =     pivot_wider(res, id_cols=id_cols, names_from=names_from, values_from=all_of(qvalue_from),     values_fn=list),
        effectsize = pivot_wider(res, id_cols=id_cols, names_from=names_from, values_from=all_of(effectsize_from), values_fn=list),
        method_order = sort(unique(res[,method_from]))
    )
}