 ###############################################################################
# Continuous value analysis

mda.continuous <- function(mda.D, continuous.cols=NULL, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(dplyr)
        require(tibble)
        require(lmerTest)})
    
    if(is.null(continuous.cols)){
        message("[MDA] mda.continuous: Undefined 'continuous.cols' parameter. You should specify this in order to use the function.")
        exit(0)
    }
    
    lmcont <- function(D, formula){

        f <- update(formula, mda.cont.col ~ .)
        
        res <- lapply(continuous.cols, function(c){
            meta_data <- D$meta_data
            meta_data[!is.na(meta_data[,c]),]
            meta_data$mda.cont.col <- meta_data[,c]
            fit <- lm(f, data=meta_data, na.action = 'na.exclude')
            s <- as.data.frame(coefficients(summary(fit)))
            s$taxa <- rep(c, dim(s)[1])
            s <- s %>% rownames_to_column("variable")
            s
        })
        res <- bind_rows(res)

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="t value"] <- "stat"
        names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"

        res
    }
    
    lmercont <- function(D, formula){
        suppressPackageStartupMessages(library(lmerTest))

        f <- update(formula, mda.group.col ~ .)

        res <- lapply(continuous.cols, function(c){
            meta_data <- D$meta_data
            meta_data[!is.na(meta_data[,c]),]
            meta_data$mda.cont.col <- meta_data[,c]
            fit <- lmer(f, data=meta_data, na.action = 'na.exclude')
            s <- as.data.frame(coefficients(summary(fit)))
            s$taxa <- rep(c, dim(s)[1])
            s <- s %>% rownames_to_column("variable")
            s
        })
        res <- bind_rows(res)

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="t value"] <- "stat"
        names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"

        res
    }

    do <- function(f_idx){
        
        fdata <- D$formula[[f_idx]]

        method <- if ( formula.ismixed(fdata$fn) ){
            lmercont
        } else { lmcont }
        
        f <- fdata$fn
        
        first_var <- formula.parts(fdata$fn.orig)[1]

        res.full <- mda.cache_load_or_run_save(D, "continuous", f_idx, method(D$count_data, fdata$data, f))

        res.full$formula <- rep(mda.deparse(fdata$fn.orig), dim(res.full)[1])
        res.full$method <- rep("lmclr", dim(res.full)[1])
        res.full <- left_join(res.full, fdata$nfreq, by="variable.mda")

        # Select only the relevant variable
        res <- res.full[startsWith(res.full$variable, mainvar),]

        res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
        res$variable <- rep(mainvar, dim(res)[1])

        return(list(res=res, res.full=res.full))
    }

    R <- lapply(1:length(D$formula$main_var), do)


    res <- bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    res.full <- bind_rows(lapply(R, function(x){x$res.full}))

    ###############################################################################
    # Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}