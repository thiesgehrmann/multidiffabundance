###############################################################################
# LMCLR

mda.lmclr <- function(mda.D, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(dplyr)
        require(tibble)
        library(lmerTest)})
    
    clr_data <- as.data.frame(scale(mda.clr(mda.relative_abundance(mda.pseudocount(D$count_data)))))

    lmclr <- function(f_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn

        method <- if ( formula.ismixed(f) ){
            lmerTest::lmer
        } else { lm }
        
        taxa <- D$nonrare

        f <- update(f, clrtaxa ~ .)
        
        data <- data.frame(fdata$data)

        res <- lapply(taxa, function(t){
            data$clrtaxa <- clr_data[rownames(data),t]
            r <- mda.trycatchempty(D, f_idx, method(f, data=data, na.action = 'na.exclude'), taxa=t)
            
            if (r$error){
                return(r$response)
            }
            fit <- r$response
            
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
        res <- dplyr::bind_rows(res)

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="t value"] <- "stat"
        names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"

        res
    }

    do <- function(f_idx){
        res.full <- mda.cache_load_or_run_save(D, f_idx, "lmclr", {lmclr(f_idx)})
        mda.common_do(D, f_idx, res.full, "lmclr")
    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)
}