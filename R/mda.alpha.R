mda.alpha <- function(mda.D, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(dplyr)
        require(tibble)
        require(vegan)
        require(lmerTest)})
    
    alpha_measure <- scale(diversity(D$count_data))

    alpha <- function(f_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn

        method <- if ( formula.ismixed(f) ){
            lmer
        } else { lm }
        
        meta_data <- data.frame(fdata$data)
        meta_data$mda.alpha <- alpha_measure

        f <- update(f, mda.alpha ~ .)
        
        r <- mda.trycatchempty(D, f_idx, method(f, data=meta_data, na.action = 'na.exclude'), taxa="mda.alpha")

        if (r$error){
            return(r$response)
        }
        fit <- r$response

        s <- as.data.frame(coefficients(summary(fit)))
        
        if (mda.isSingular(fit)){
            s[,"Pr(>|t|)"] <- NA
            s[,"comment"] <- paste0(c("Rank deficient: singular.", r$message), collapse='\n')
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
        res.full <- mda.cache_load_or_run_save(D, f_idx, "alpha", {alpha(f_idx)})
        mda.common_do(D, f_idx, res.full, "alpha", skip_taxa_sel=TRUE)
    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)

}