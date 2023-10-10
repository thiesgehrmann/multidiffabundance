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
    
    group <- function(f_idx){

        fdata <- D$formula[[f_idx]]
        f <- fdata$fn

        method <- if ( formula.ismixed(f) ){
            glmer
        } else { glm }
        
        f <- update(f, mda.group.value ~ .)
        
        all.res <- lapply(group.cols, function(col){
            rel_idx <-    !is.na(D$meta_data[,col])
            rel_meta <-   data.frame(fdata$data[rel_idx, ])
            rel_meta$mda.group.col <- D$meta_data[rel_idx, col]
            
            indiv.res <-  lapply(unique(rel_meta$mda.group.col), function(g){
                
                indiv.meta_data <- rel_meta
                indiv.meta_data$mda.group.value <- as.numeric(rel_meta$mda.group.col == g)
                
                r <- mda.trycatchempty(D, f_idx, method(f, data=indiv.meta_data, na.action = 'na.exclude', family="binomial"), taxa=g)
            
                if (r$error){
                    return(r$response)
                }
                fit <- r$response

                
                s <- as.data.frame(coefficients(summary(fit)))
                if (mda.isSingular(fit)){
                    s[,"Pr(>|z|)"] <- NA
                    s[,"comment"] <- paste0(c("Rank deficient: singular.", r$message), collapse='\n')
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
        res.full <- mda.cache_load_or_run_save(D, f_idx, "group", {group(f_idx)})
        mda.common_do(D, f_idx, res.full, "group", skip_taxa_sel=TRUE)
    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)
}