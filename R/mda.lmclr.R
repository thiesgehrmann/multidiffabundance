 ###############################################################################
# LMCLR

mda.lmclr <- function(mda.D, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(dplyr)
        require(tibble)
        library(lme4)})
    
    clr_data <- as.data.frame(scale(mda.clr(mda.relative_abundance(mda.pseudocount(D$count_data)))))

    lmclr <- function(count_data, meta_data, formula, taxa=NULL, method){
        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        f <- update(formula, clrtaxa ~ .)

        res <- lapply(taxa, function(t){
            meta_data$clrtaxa <- clr_data[,t]
            r <- tryCatch({
                    list(fit=method(f, data=meta_data, na.action = 'na.exclude'), error=FALSE)
                },
                error=function(err){
                    return(list(fit=mda.empty_output(fdata, err$message), error=TRUE))
                })

            if (r$error){
                return(r$fit)
            }
            fit <- r$fit
            
            s <- as.data.frame(coefficients(summary(fit)))
            s[,"comment"] <- NA
            if (mda.isSingular(fit)){
                s[,"Pr(>|t|)"] <- NA
                s[,"comment"] <- "Rank deficient: singular"
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
        f <- fdata$fn

        method <- if ( formula.ismixed(f) ){
            lmer
        } else { lm }

        res.full <- mda.cache_load_or_run_save(D, "lmclr", f_idx, lmclr(D$count_data, fdata$data, f, D$nonrare, method=method))

        mda.common_do(D, res.full, "lmclr", fdata)
    }

    R <- lapply(1:length(D$formula), do)

    mda.common_output(R)
}
