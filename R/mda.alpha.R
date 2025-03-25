#' @importFrom lmerTest lmer
#' @importFrom stringr str_replace_all
#' @export
mda.alpha <- function(mda.D, alpha.index=c("shannon"), ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(dplyr)
        require(tibble)
        require(vegan)})
    
    index <- intersect(alpha.index, c("shannon", "simpson", "invsimpson"))
    
    if(length(index) == 0){
        message("[MDA] mda.alpha: No valid index value(s) defined. Options must be [shannon, simpson, invsimpson].")
        exit(0)
    }
    
    alpha_indices <- as.data.frame(lapply(index, function(i){scale(diversity(D$count_data, index=i))}))
    colnames(alpha_indices) <- index

    alpha <- function(f_idx, i_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn

        method <- if ( formula.ismixed(f) ){
            lmer
        } else { lm }
        
        meta_data <- data.frame(fdata$data)
        meta_data$mda.alpha <- alpha_indices[rownames(meta_data),index[i_idx]]

        taxa_name <- paste0(c("mda","alpha",index[i_idx]), collapse='.')

        f <- update(f, mda.alpha ~ .)

        r <- mda.trycatchempty(D, f_idx, method(f, data=meta_data, na.action = 'na.exclude'), taxa=taxa_name)

        if (r$error){
            return(r$response)
        }
        fit <- r$response

        s <- as.data.frame(coefficients(summary(fit)))
        s$comment <- ""

        if (mda.isSingular(fit)){
            s[,"Pr(>|t|)"] <- NA
            s[,"comment"] <- str_replace_all(paste0(c("Rank deficient: singular.", r$message), collapse=' ; '),'\n',' ; ')
        }
        
        s$taxa <- taxa_name
        s <- s %>% rownames_to_column("variable.mda")
        res <- s

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="t value"] <- "stat"
        names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"

        res
    }

    do <- function(f_idx){
        res.full <- lapply(1:length(index), function(i_idx){
            m <- paste(c("alpha",index[i_idx]), collapse=".")
            mda.cache_load_or_run_save(D, f_idx, m, {alpha(f_idx, i_idx)})
            })
        res.full <- bind_rows(res.full)
        mda.common_do(D, f_idx, res.full, "alpha", skip_taxa_sel=TRUE)

    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)

}
