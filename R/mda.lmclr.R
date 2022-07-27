 ###############################################################################
# LMCLR

mda.lmclr <- function(mda.D){
    D <- mda.D
    
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tibble))
    
    clr_data <- as.data.frame(scale(mda.clr(mda.relative_abundance(mda.pseudocount(D$count_data)))))

    lmclr <- function(count_data, meta_data, formula, mainvar, taxa=NULL){
        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        f <- update(formula, clrtaxa ~ .)

        res <- lapply(taxa, function(t){
            meta_data$clrtaxa <- clr_data[,t]
            fit <- lm(f, data=meta_data, na.action = 'na.exclude')
            s <- as.data.frame(coefficients(summary(fit)))
            s$taxa <- rep(t, dim(s)[1])
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
    
    lmerclr <- function(count_data, meta_data, formula, mainvar, taxa=NULL){
        suppressPackageStartupMessages(library(lmerTest))
        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        f <- update(formula, clrtaxa ~ .)

        res <- lapply(taxa, function(t){
            meta_data$clrtaxa <- clr_data[,t]
            fit <- lmer(f, data=meta_data, na.action = 'na.exclude')
            s <- as.data.frame(coefficients(summary(fit)))
            s$taxa <- rep(t, dim(s)[1])
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

        method <- if ( (length(D$formula$rand_intercept[[f_idx]]) + length(D$formula$rand_slope[[f_idx]])) > 0 ){
            lmerclr
        } else { lmclr }
        f <- D$formula$formula[[f_idx]]
        mainvar <- D$formula$main_var[f_idx]

        res.full <- mda.cache_load_or_run_save(D, "lmclr", f, method(D$count_data, D$meta_data, f, mainvar, D$nonrare))

        res.full$formula <- rep(mda.deparse(f), dim(res.full)[1])
        res.full$method <- rep("lmclr", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

        # Select only the relevant variable ( taxa are selected already in lmclr, but we repeat it here for safety )
        res <- res.full[res.full$taxa %in% D$nonrare,]
        res <- res[startsWith(res$variable, mainvar),]

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