mda.group <- function(mda.D, group.cols=NULL, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(dplyr)
        require(tibble)
        require(lmerTest)})
    
    if(is.null(group.cols)){
        message("[MDA] mda.group: Undefined 'group.cols' parameter. You should specify this in order to use the function.")
        exit(0)
    }
    
    group <- function(count_data, meta_data, formula, method){

        f <- update(formula, mda.group.value ~ .)
        
        all.res <- lapply(group.cols, function(col){
            rel_idx <-    !is.na(D$meta_data[,col])
            rel_meta <-   meta_data[rel_idx, ]
            rel_meta$mda.group.col <- D$meta_data[rel_idx, col]
            
            indiv.res <-  lapply(unique(rel_meta$mda.group.col), function(g){
                
                indiv.meta_data <- rel_meta
                meta_data$mda.group.value <- as.numeric(rel_meta$mda.group.col == g)

                r <- tryCatch({
                        list(fit=method(f, data=meta_data, na.action = 'na.exclude', family="binomial"), error=FALSE)
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
                    s[,"Pr(>|z|)"] <- NA
                }
                s$taxa <- rep(paste0(c(col,".", g), collapse=""), dim(s)[1])
                s <- s %>% rownames_to_column("variable.mda")
                s
            })
            indiv.res <- bind_rows(indiv.res)
        })
        
        res <- bind_rows(all.res)

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="z value"] <- "stat"
        names(res)[names(res)=="Pr(>|z|)"] <- "pvalue"

        res
    }


    do <- function(f_idx){
        
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn

        method <- if ( formula.ismixed(f) ){
            glmer
        } else { glm }

        res.full <- mda.cache_load_or_run_save(D, "group", f_idx, group(D$count_data, fdata$data, f, method=method))

        mda.common_do(D, res.full, "group", fdata, skip_taxa_sel=TRUE)

    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)

}
