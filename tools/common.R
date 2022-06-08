library("tidyverse")

mda.load <- function(args){
    
    if (length(args) != 4){
        stop("Too few arguments. There should be 4!")
    }
    
    abundance <- args[1]
    meta <- args[2]
    formula.data <- args[3]
    outprefix <- args[4]

    ###############################################################################
    # Prepare formula

    F <- c()

    raw <- if (file.exists(formula.data)){
        unlist(strsplit(read_file(formula.data), "\n"))
    } else {
        formula.data
    }

    form <- lapply(raw, as.formula)

    fterms <- lapply(form, function(x){labels(terms(x))})
    fterms_fixed  <- lapply(fterms, function(x){ x[!grepl("\\|", x)]})
    fterms_random <- lapply(fterms, function(x){ x[grepl("\\|", x)]})

    F$raw <- raw
    F$formula <- form
    F$main_var <- unlist(lapply(fterms_fixed, function(x){unlist(x[1])}))
    F$adj_vars <- unlist(lapply(fterms_fixed, function(x){paste0(x[-1], collapse="+")}))
    F$rand_vars <- unlist(lapply(fterms_random, function(ftr){
        paste0(lapply(ftr, function(x){paste0(c('(', x, ')'), collapse="")}), collapse="+")
    }))
    F$norand <- lapply(fterms_fixed, function(ftf){
        as.formula(paste0(c("~",paste0(ftf, collapse="+")),collapse=""))
    })

    if (str_length(paste0(unlist(fterms_random), collapse="")) > 0){
        print("No mixed effect terms are allowed in this implementation. Random effects will be ignored. Please format your random effects as fixed effects. Sorry.")
    }

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
                                                        
    ###############################################################################
    
    dat <- c()
    dat$count_data  <- count_data
    dat$nonrare     <- nonrare
    dat$meta_data   <- meta_data
    dat$formula     <- F
    dat$outprefix   <- outprefix
    
    return(dat)
    
}

###############################################################################

mda.nonrare_taxa <- function(table , cutoff_pct) {
    cutoff  <- ceiling(cutoff_pct * nrow(table))
    nonzero <- colnames(table)[colSums(table > 0) >= cutoff]
    return(nonzero)
}

###############################################################################

mda.relative_abundance <- function(count_data){
    depth <- rowSums(count_data)
    radata <- apply(count_data, 2, function(x){x/depth})
    return(radata)
}

###############################################################################

mda.pseudocount <- function(count_data){
    pcount <- (rowSums(count_data) / max(rowSums(count_data)))
    pdata <- apply(count_data,2,function(x){x+pcount})
    return(pdata)
}

###############################################################################
                                                      
mda.clr <- function(df){
    denom <- exp(rowMeans(log(df)))
    df / denom
}

###############################################################################

mda.meta.n <- function(D, var){
    sum(unlist(lapply(D$meta_data[,var],function(x){!(is.na(x))})))
}

###############################################################################

mda.meta.freq <- function(D, var){
    if(is.numeric(D$meta_data[,var])){
        ""
    } else {
        paste0(mapply(function(x,y){paste0(c(x,y), collapse=":")}, names(table(D$meta_data[,var])), as.character(table(D$meta_data[,var]))), collapse=", ")
    }
}