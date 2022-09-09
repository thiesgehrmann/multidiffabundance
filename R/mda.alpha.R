###############################################################################
# ALPHA

mda.alpha <- function(mda.D, ...){
    D <- mda.D
    
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tibble))
    suppressPackageStartupMessages(require(vegan))
    
    D$meta_data$mda.alpha <- scale(diversity(D$count_data))

    lmalpha <- function(count_data, meta_data, formula, mainvar, taxa=NULL){
        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        f <- update(formula, mda.alpha ~ .)

        fit <- lm(f, data=meta_data, na.action = 'na.exclude')
        s <- as.data.frame(coefficients(summary(fit)))
        s$taxa <- c("mda.alpha")
        s <- s %>% rownames_to_column("variable")
        res <- s

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="t value"] <- "stat"
        names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"

        res
    }
    
    lmeralpha <- function(count_data, meta_data, formula, mainvar, taxa=NULL){
        suppressPackageStartupMessages(require(lmerTest))
        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        f <- update(formula, mda.alpha ~ .)
        fit <- lmer(f, data=meta_data, na.action = 'na.exclude')
        s <- as.data.frame(coefficients(summary(fit)))
        s$taxa <- c("mda.alpha")
        s <- s %>% rownames_to_column("variable")
        res <- s

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="t value"] <- "stat"
        names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"

        res
    }

    do <- function(f_idx){

        method <- if ( (length(D$formula$rand_intercept[[f_idx]]) + length(D$formula$rand_slope[[f_idx]])) > 0 ){
            lmeralpha
        } else { lmalpha }
        f <- D$formula$formula[[f_idx]]
        mainvar <- D$formula$main_var[f_idx]

        res.full <- mda.cache_load_or_run_save(D, "alpha", f, method(D$count_data, D$meta_data, f, mainvar, D$nonrare))

        res.full$formula <- rep(mda.deparse(f), dim(res.full)[1])
        res.full$method <- rep("alpha", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

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