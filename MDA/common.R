###############################################################################
# MDA Common R functions
###############################################################################


###############################################################################
# Input processing function

mda.load <- function(args){
    require(tools)
    require("tidyverse")
    
    if (length(args) != 4){
        stop("Too few arguments. There should be 4!")
    }
    
    abundance <- args[1]
    meta <- args[2]
    formula.data <- args[3]
    outprefix <- args[4]

    ###############################################################################
    # Prepare formula

    raw.formula.data <- mda.load_formula_input(formula.data)
    if (!mda.verify_formula_input(raw.formula.data)){
        print("Error: Errors in processing formula input.")
        quit(1)
    }
    FD <- mda.process_formula_input(raw.formula.data)

    ###############################################################################
    # Load the data

    count_data <- as.data.frame(read_tsv(abundance))
    row.names(count_data) <- count_data[,1]
    count_data <- count_data[,-1]

    nonrare <- mda.nonrare_taxa(count_data, 0.1) 

    meta_data <- as.data.frame(read_tsv(meta))
    row.names(meta_data) <- meta_data[,1]
    meta_data <- meta_data[,-1]
    numeric_meta <- colnames(meta_data)[unlist(lapply(colnames(meta_data), function(x)is.numeric(meta_data[,x])))]
    meta_data[,numeric_meta] <- scale(meta_data[,numeric_meta])
                                                      
    meta_data <- meta_data[rownames(count_data), ]
                                                        
    ###############################################################################
    
    dat <- c()
    dat$count_data  <- count_data
    dat$nonrare     <- nonrare
    dat$meta_data   <- meta_data
    dat$formula     <- FD
    dat$outprefix   <- outprefix
    dat$cacheprefix <- paste0(c(outprefix, "mda.cache", paste0(md5sum(c(abundance, meta)), collapse=".")), collapse="/")
    print(dat$cacheprefix)
    mda.mkdirp(dirname(dat$cacheprefix))
    
    return(dat)
    
}
                                                      
mda.mkdirp <- function(dir){
    if (!dir.exists(dir)){
        mda.mkdirp(dirname(dir))
        dir.create(dir)
    }
}

###############################################################################
# Formula processing functions
                                                      
mda.load_formula_input <- function(formula_input){
  raw <- if (file.exists(formula_input)){
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
  
  fterms <- lapply(form, function(x){labels(terms(x))})
  fterms_fixed  <- lapply(fterms, function(x){ x[!grepl("\\|", x)]})
  fterms_random <- lapply(fterms, function(x){ x[grepl("\\|", x)]})
  
  FD <- c()
  FD$raw <- raw_formula
  FD$formula <- form
  FD$main_var <- unlist(lapply(fterms_fixed, function(x){unlist(x[1])}))
  FD$adj_vars <- unlist(lapply(fterms_fixed, function(x){paste0(x[-1], collapse="+")}))
  FD$rand_vars <- unlist(lapply(fterms_random, function(ftr){
    paste0(lapply(ftr, function(x){paste0(c('(', x, ')'), collapse="")}), collapse="+")
  }))
  FD$norand <- lapply(fterms_fixed, function(ftf){
    as.formula(paste0(c("~",paste0(ftf, collapse="+")),collapse=""))
  })
  
  if (str_length(paste0(unlist(fterms_random), collapse="")) > 0){
    print("WARNING: No mixed effect terms are allowed in this implementation. Random effects will be ignored. Please format your random effects as fixed effects. Sorry.")
  }
  
  FD
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
  form.fmt <- tolower(format(form))
  form.fmt <- gsub("[~():._! ]", "", form.fmt)
  form.fmt <- gsub("[*]", "M", form.fmt)
  form.fmt <- gsub("[+]", "P", form.fmt)
  form.fmt <- gsub("[:]", "I", form.fmt)
  form.fmt <- gsub("[|]", "R", form.fmt)
  form.fmt <- gsub("[-]", "S", form.fmt)
  form.fmt
  
  filename <- paste0(c(paste0(c(outprefix, method, form.fmt), collapse=collapse), suffix), collapse=".")
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
                                                      
mda.cache_load_or_run_save <- function(outputprefix, method, form, expr) {
    cache.file <- mda.cache_filename(outputprefix, method, form, suffix="rds")
    data <- if (file.exists(cache.file)){
        print(paste0(c("Loading:", cache.file, "for", method, ",", format(form)), collapse=" "))
        readRDS(cache.file)
    } else{
        data <- expr
        print(paste0(c("Executing & storing:", cache.file, "for", method, ",", format(form)), collapse=" "))
        saveRDS(data, cache.file)
        data
    }
    data
}
                                                      
###############################################################################