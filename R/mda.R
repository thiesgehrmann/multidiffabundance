###############################################################################
# MDA Common R functions
###############################################################################
                                                      
###############################################################################
# Input loading functions
#' @export
mda.from_cmdargs <- function(args, ...){
    suppressPackageStartupMessages({
        require(tools)
        require("tidyverse")})
    
    if (length(args) != 4){
        mda.message("Too few arguments. There should be 4!", type="error")
        quit(1)
    }
    
    abundance <- args[1]
    meta <- args[2]
    formula.data <- args[3]
    outprefix <- args[4]

    mda.from_files(abundance, meta, formula.data, outprefix, ...)
}
        
#' @export                                              
mda.from_files <- function(abundance, meta, formula.data, outprefix=tempdir(), ...){
    # Prepare formula
    suppressPackageStartupMessages(require(tidyverse))

    raw.formula.data <- mda.load_formula_input(formula.data)
    if (!mda.verify_formula_input(raw.formula.data)){
        mda.message("Errors in processing formula input.", type="error")
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
#' Create an mda object from a tidytacos object and a formula
#' @export
mda.from_tidytacos <- function(ta, formulas, output_dir=tempdir(), ...){
    suppressPackageStartupMessages({
        require(dplyr)
        require(tidyr)})
    meta_data <- as.data.frame(ta$samples)
    rownames(meta_data) <- meta_data$sample_id
    
    count_data <- as.data.frame(pivot_wider(ta$counts, id_cols="sample_id", names_from="taxon_id", values_from="count", values_fill=0))
    rownames(count_data) <- count_data$sample_id
    count_data <- count_data[,-1]

    union.samples = sort(union(rownames(count_data), rownames(meta_data)))
    intersect.samples = sort(intersect(rownames(count_data), rownames(meta_data)))
    
    if ( length(union.samples) != length(intersect.samples) ){
        mda.message("tidytacos object has differing samples present in abundance and meta data!")
    }
    
    if (!"taxon" %in% colnames(ta$taxa)) {
        ta <- ta %>% tidytacos::add_taxon_name()
        ta$taxa <- ta$taxa %>%
          dplyr::rename(taxon=taxon_name) #needed by ancombc2
    }
    
    dat <- mda.create(count_data[intersect.samples,], meta_data[intersect.samples,], formulas, output_dir, ...)
    dat$ta <- ta
    
    return(dat)
}

#' @export
mda.from_phyloseq <- function(phys, formulas, output_dir=tempdir(), ...){
    mda.message("UNIMPLEMENTED", type="error")
}
#' @export                                            
mda.create <- function(count_data, meta_data, formulas, output_dir=tempdir(), usecache=TRUE, recache=FALSE, nonrare.pct=0.1, ...){
    suppressPackageStartupMessages(require(digest))
    
    if (! all(rownames(count_data) == rownames(meta_data))) {
        mda.message("Rownames of meta data and count data do not match!", type="error")
        return(NULL)
    }
    
    nonrare <- mda.nonrare_taxa(count_data, nonrare.pct) 

    numeric_meta <- colnames(meta_data)[unlist(lapply(colnames(meta_data), function(x)is.numeric(meta_data[,x])))]
    meta_data[,numeric_meta] <- scale(meta_data[,numeric_meta])

    FD <- lapply(as.list(formulas), function(f){formula.reformulate(as.formula(f), meta_data)})
                                                      
    
    dat <- c()
    dat$count_data  <- count_data
    dat$nonrare     <- nonrare
    dat$nonrare.pct <- nonrare.pct
    dat$meta_data   <- meta_data
    dat$formula     <- FD
    dat$outprefix   <- output_dir
    dat$usecache    <- usecache
    dat$recache     <- recache
    checksum <- digest(count_data, algo="md5")
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
                                                      
mda.as_tidytacos <- function(D){
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

#' @export                                                      
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
#' @export
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

#' Permutes over all possible orders of effects in a given formula and returns a list of the possibilities.
#'
#' @param form a formula to permute over, eg: ~ a:b + a + b.
#'
#' @export
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
# output cache functions

mda.cache_load_or_run_save <- function(mda.D, f_idx, method, expr, order_invariant=TRUE, extra=NULL) {
    D <- mda.D
    
    checksum <- if (is.null(f_idx)){ "" } else { if (order_invariant){ D$formula[[f_idx]]$checksum.order_invariant  } else {  D$formula[[f_idx]]$checksum.order_variant } }
    mainvar  <- if (is.null(f_idx)){ "mda_crossmethod_exec" } else { formula.parts(D$formula[[f_idx]]$fn.orig)[1] }

    checksum <- if(is.null(extra)){checksum} else { digest(c(checksum, extra), algo="md5") }
    
    cache.file <- paste0(c(paste0(c(D$cacheprefix, method, checksum), collapse='_'), "rds"), collapse=".")
    
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
                  
###############################################################################
# trycatch functions
                  
mda.trycatchempty <- function(D, f_idx, expr, taxa=NA, empty=NULL){
    fdata <- D$formula[[f_idx]]
    
    data <- tryCatch(withCallingHandlers({
            message <- NULL
            list(response=expr, error=FALSE, message=message)
          }, 
          warning = function(warn) {
            message <<- warn$message
            invokeRestart("muffleWarning")
          }
        ),
        error=function(err){
            response <- if (is.null(empty)) { mda.empty_output(D, f_idx, comment=err$message, taxa=taxa) } else { empty }
            return(list(response=response, error=TRUE, message=err$message))
        }
    )
    data
}

###############################################################################
# Empty output format
                  
mda.empty_output <- function(D, f_idx, comment=NA, taxa=NA){
    fdata <- D$formula[[f_idx]]
    empty.res <- fdata$nfreq[fdata$parts.fixed, c('variable.mda'),drop=FALSE]

    empty.res[,c('se','taxa','pvalue','effectsize','df','stat')] <- NA
    
    empty.res <- if ((length(taxa) == 0) | all(is.na(taxa)) ) {
        empty.res
    } else {
        bind_rows(lapply(taxa, function(t){
            e <- data.frame(empty.res)
            e$taxa <- t
            e
        } ))
    }

    empty.res <- transform(empty.res,
                           se = as.numeric(se),
                           taxa = as.character(taxa),
                           pvalue = as.numeric(pvalue),
                           effectsize = as.numeric(effectsize),
                           df = as.numeric(df),
                           stat = as.numeric(stat))
    empty.res$comment <- comment
    empty.res
}
###############################################################################

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
#' Summary function
#' @export
mda.summary <- function(res, id_cols = "taxa", names_from = "variable", method_from = "method",
    pvalue_from = "pvalue", qvalue_from = "qvalue", effectsize_from = "effectsize", 
    values_fn = list, pvalue_threshold = 0.05, qvalue_threshold = 0.05,
    method.groups = list(
        "da" = c("aldex2", "ancombc2",'corncob','deseq2','corncob','limma','lmclr','maaslin2','zicoseq'),
        "alpha" = c("alpha"),
        "beta"  = c("beta"),
        "group" = c("group"),
        "continuous" = c("continuous"))) {
    
    res <- as.data.frame(res)

    mda.summary.internal <- function (res)
    {
        res <- dplyr::arrange(res, method_from, id_cols, names_from)
        methods <- sort(unique(res[,method_from]))

        list(
            nsig =       tidyr::pivot_wider(res, id_cols=all_of(id_cols), names_from=all_of(names_from), values_from=all_of(qvalue_from),     values_fn=function(v){as.integer(sum(v < qvalue_threshold))}),
            pvalue =     lapply(methods, function(method) { tidyr::pivot_wider(res[res[,method_from] == method,], id_cols=all_of(id_cols), names_from=all_of(names_from), values_from=all_of(pvalue_from),     values_fn=values_fn)}),
            qvalue =     lapply(methods, function(method) { tidyr::pivot_wider(res[res[,method_from] == method,], id_cols=all_of(id_cols), names_from=all_of(names_from), values_from=all_of(qvalue_from),     values_fn=values_fn)}),
            effectsize = lapply(methods, function(method) { tidyr::pivot_wider(res[res[,method_from] == method,], id_cols=all_of(id_cols), names_from=all_of(names_from), values_from=all_of(effectsize_from), values_fn=values_fn)}),
            methods    = methods
        )
    }
    
    group.summary <- lapply(method.groups, function(group){
        group.res <- res[res[,method_from] %in% group,]
        if(nrow(group.res) > 0) {mda.summary.internal(group.res)} else {NULL}
    })

    group.summary
    
}
                  
###############################################################################
#' A function to merge the outputs of several mda runs
#' @export
mda.merge_results <- function(res.list){

    comb <- list()

    comb$res      <- bind_rows(lapply(res.list, function(x){x$res}))
    comb$res.full <- bind_rows(lapply(res.list, function(x){x$res.full}))
    comb$summary  <- mda.summary(comb$res)
    
    comb
}


###############################################################################
# A wrapper to be able to check singularity for both lmer and lm objects

mda.isSingular <- function(fit){
    suppressPackageStartupMessages({
        require(lmerTest)})

    if (length(intersect(class(fit), c("lmerModLmerTest", "glmerMod"))) > 0){
        isSingular(fit)
    }
    else if (length(intersect(class(fit), c("lm", "glm"))) > 0){
        p <- length(attributes(fit$terms)$term.labels) + attributes(fit$terms)$intercept
        fit$rank < p
    }
    else { FALSE } # We don't know what to do in this case...
}

###############################################################################
# Common post-do formatting
#' @importFrom stringr str_replace
#' @importFrom dplyr left_join
mda.common_do <- function(D, f_idx, res.full, method, skip_taxa_sel=FALSE){
    print(1)
    fdata <- D$formula[[f_idx]]
    
    res.full$formula <- rep(mda.deparse(fdata$fn.orig), dim(res.full)[1])
    res.full$method <- rep(method, dim(res.full)[1])
    res.full <- left_join(res.full, fdata$nfreq, by="variable.mda")

    # Select only the relevant taxa ( taxa are selected already in some methods, but we repeat it here for safety )
    res.full <- if(skip_taxa_sel){
        res.full 
    } else {
        missing <- if (length(setdiff(D$nonrare, res.full$taxa)) == 0){
                NULL
            } else {
                m <- mda.empty_output(D, f_idx, comment="This method did not test this taxa.", taxa=setdiff(D$nonrare, res.full$taxa))
                m$formula <- rep(mda.deparse(fdata$fn.orig), dim(m)[1])
                m$method <- method
                left_join(m, fdata$nfreq, by="variable.mda")
            }
        
        res.full <- res.full[res.full$taxa %in% D$nonrare,]
        if (is.null(missing)){ res.full } else { bind_rows(res.full, missing) }
    }
    
    if (is.null(res.full$comment)){
        res.full$comment <- ""
    } else {
        res.full$comment <- lapply(res.full$comment, function(x){str_replace(x, '\n', '')})
    }

    # Select only the first variable in the original function
    first_var <- formula.parts(fdata$fn.orig)[1]
    res <- res.full[res.full$variable == first_var,]
    res <- res[!is.na(res$variable),]
    
    # FDR correction
    res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
    res.full$qvalue.withinformula <- p.adjust(res.full$pvalue, "fdr")

    return(list(res=res, res.full=res.full))
}
                  
###############################################################################
# Common output format
                  
mda.common_output <- function(R){

    column.initial <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","formula","method","n","freq","comment")
    column.order   <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq","comment")

    
    res <- dplyr::bind_rows(lapply(R, function(x){
        v <- x$res[,intersect(column.initial, colnames(x$res))]
        v[,setdiff(column.initial, colnames(x$res))] <- NA
        v
    }))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    rownames(res) <- NULL
    
    res.full <- dplyr::bind_rows(lapply(R, function(x){
        v <- x$res.full[,intersect(column.initial, colnames(x$res.full))]
        v[,setdiff(column.initial, colnames(x$res.full))] <- NA
        v
    }))
    res.full$qvalue <- p.adjust(res.full$pvalue, "fdr")
    rownames(res.full) <- NULL

    # Select output columns

    res <-          res[,column.order]
    res.full <-     res.full[,column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}
###############################################################################
# Error/Warning display in terminal

mda.message <- function(msg, type="warning") {
    if (type != "warning") { type <- "error" }
    type <- toupper(type)
    func <- rlang::caller_call()[[1]]
    msg <- paste0(c("[MDA] ", func, " ", type, ": ", msg))
    message(msg)
}                  

