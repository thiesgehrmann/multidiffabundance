 ###############################################################################
# LMCLR

mda.lmclr <- function(mda.D, ...){
    D <- mda.D
    
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tibble))
    library(lme4)
    
    clr_data <- as.data.frame(scale(mda.clr(mda.relative_abundance(mda.pseudocount(D$count_data)))))

    lmclr <- function(count_data, meta_data, formula, taxa=NULL){
        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        f <- update(formula, clrtaxa ~ .)

        res <- lapply(taxa, function(t){
            meta_data$clrtaxa <- clr_data[,t]
            fit <- lm(f, data=meta_data, na.action = 'na.exclude')
            
            s <- as.data.frame(coefficients(summary(fit)))
            #if (isSingular(fit)){
            #    s[,"Pr(>|t|)"] <- NA
            #    }
            
            s$taxa <- rep(t, dim(s)[1])
            s <- s %>% rownames_to_column("variable.mda")
            s
        })
        res <- bind_rows(res)

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="t value"] <- "stat"
        names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"

        res
    }
    
    lmerclr <- function(count_data, meta_data, formula, taxa=NULL){
        suppressPackageStartupMessages(library(lmerTest))
        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        f <- update(formula, clrtaxa ~ .)

        res <- lapply(taxa, function(t){
            meta_data$clrtaxa <- clr_data[,t]
            fit <- lmer(f, data=meta_data, na.action = 'na.exclude')
            s <- as.data.frame(coefficients(summary(fit)))
            if (isSingular(fit)){
                s[,"Pr(>|t|)"] <- NA
                }
            s$taxa <- rep(t, dim(s)[1])
            s <- s %>% rownames_to_column("variable.mda")
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
            lmerclr
        } else { lmclr }
        
        f <- fdata$fn
        
        first_var <- formula.parts(fdata$fn.orig)[1]

        res.full <- mda.cache_load_or_run_save(D, "lmclr", f_idx, method(D$count_data, fdata$data, f, D$nonrare))

        res.full$formula <- rep(mda.deparse(fdata$fn.orig), dim(res.full)[1])
        res.full$method <- rep("lmclr", dim(res.full)[1])
        res.full <- left_join(res.full, fdata$nfreq, by="variable.mda")

        # Select only the relevant variable ( taxa are selected already in lmclr, but we repeat it here for safety )
        res.full <- res.full[res.full$taxa %in% D$nonrare,]
        
        res <- res.full[res.full$variable == first_var,]
        res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
        
        res.full$qvalue.withinformula <- p.adjust(res.full$pvalue, "fdr")

        return(list(res=res, res.full=res.full))
    }

    R <- lapply(1:length(D$formula), do)


    res <- bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    
    res.full <- bind_rows(lapply(R, function(x){x$res.full}))
    res.full$qvalue <- p.adjust(res.full$pvalue, "fdr")

    ###############################################################################
    # Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    res.full <- res.full[, column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}