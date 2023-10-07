 ###############################################################################
# Continuous value analysis

mda.continuous <- function(mda.D, continuous.cols=NULL, continuous.scale=TRUE, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(dplyr)
        require(tibble)
        require(lmerTest)})
    
    if(is.null(continuous.cols)){
        message("[MDA] mda.continuous: Undefined 'continuous.cols' parameter. You should specify this in order to use the function.")
        exit(0)
    }
    
    continuous <- function(count_data, meta_data, formula, method){

        f <- update(formula, mda.cont.col ~ .)
        
        res <- lapply(continuous.cols, function(c){
            meta_data$mda.cont.col <- if(continuous.scale) { as.vector(scale(D$meta_data[,c]), ) } else { D$meta_data[,c] }

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
            s$taxa <- rep(c, dim(s)[1])
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

        res.full <- mda.cache_load_or_run_save(D, "continuous", f_idx, continuous(D$count_data, fdata$data, f, method), extra=continuous.cols)

        mda.common_do(D, res.full, "continuous", fdata, skip_taxa_sel=TRUE)

    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)

}
