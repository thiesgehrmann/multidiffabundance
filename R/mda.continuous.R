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
    
    continuous <- function(f_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn

        method <- if ( formula.ismixed(f) ){
            lmer
        } else { lm }
        f <- update(f, mda.cont.col ~ .)
        
        meta_data <- data.frame(fdata$data)
        res <- lapply(continuous.cols, function(c){
            meta_data$mda.cont.col <- if(continuous.scale) { as.vector(scale(D$meta_data[,c]), ) } else { D$meta_data[,c] }
            
            r <- mda.trycatchempty(D, f_idx, method(f, data=meta_data, na.action = 'na.exclude'), taxa=c)

            if (r$error){
                return(r$response)
            }
            fit <- r$response

            s <- as.data.frame(coefficients(summary(fit)))
            if (mda.isSingular(fit)){
                s[,"Pr(>|z|)"] <- NA
                s[,"comment"] <- paste0(c("Rank deficient: singular.", r$message), collapse='\n')
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
        res.full <- mda.cache_load_or_run_save(D, f_idx, "continuous", {continuous(f_idx)}, extra=continuous.cols)
        mda.common_do(D, f_idx, res.full, "continuous", skip_taxa_sel=TRUE)
    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)

}