mda.alpha <- function(mda.D, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(dplyr)
        require(tibble)
        require(vegan)})
    
    alpha_measure <- scale(diversity(D$count_data))

    alpha <- function(count_data, meta_data, formula, mainvar, taxa=NULL, method){
        meta_data$mda.alpha <- alpha_measure
        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        f <- update(formula, mda.alpha ~ .)

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
        
        if (mda.isSingular(fit)){
                    s[,"Pr(>|t|)"] <- NA
        }
        s$taxa <- c("mda.alpha")
        s <- s %>% rownames_to_column("variable.mda")
        res <- s

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

        res.full <- mda.cache_load_or_run_save(D, "alpha", f_idx, alpha(D$count_data, fdata$data, f, D$nonrare, method=method))
        mda.common_do(D, res.full, "alpha", fdata, skip_taxa_sel=TRUE)
    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)

}