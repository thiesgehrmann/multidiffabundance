#' @importFrom stringr str_replace
#' @importFrom dplyr rename
#' @importFrom tidyr gather
#' @export
mda.ancombc2 <- function(mda.D, ...){
    suppressPackageStartupMessages({
        require("ANCOMBC")
        require(tidyr)
        require(dplyr)
        require(tibble)
        require(stringr)})
    
    D <- mda.D

    OTU <- phyloseq::otu_table(t(D$count_data), taxa_are_rows = T)

    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]

        f.fixed <- paste0(fdata$parts.fixed, collapse=' + ')
        f.rand  <- if (formula.ismixed(fdata$fn)) { paste0(lapply(fdata$parts.random, function(v){paste0(c('(',v,')'), collapse='')}), collapse='+') } else { NULL }

        sampledata <- phyloseq::sample_data(fdata$data, errorIfNULL = F)
        phylo <- phyloseq::merge_phyloseq(OTU, sampledata)

        r <- mda.trycatchempty(D, f_idx, {
            mda.cache_load_or_run_save(D, f_idx, "ancombc2", 
                    ANCOMBC::ancombc2(data = phylo, fix_formula = f.fixed, rand_formula = f.rand, 
                            p_adj_method = "holm", prv_cut=0, lib_cut = 1000, 
                            struc_zero = FALSE, global = FALSE, alpha = 0.05, ...) )
        }, taxa=D$nonrare)
            
        if (r$error){
            return(mda.common_do(D, f_idx, r$response, "ancombc2", skip_taxa_sel = FALSE))
        }

        fit = r$response$res
        coeff <- fit[,c("taxon", colnames(fit)[startsWith(colnames(fit), 'lfc_')])]
        colnames(coeff) <- sapply(colnames(coeff), function(x){str_replace(x,"lfc_","")})
        
        se <- fit[,c("taxon", colnames(fit)[startsWith(colnames(fit), 'se_')])]
        colnames(se) <- sapply(colnames(se), function(x){str_replace(x,"se_","")})
        
        stat <- fit[,c("taxon", colnames(fit)[startsWith(colnames(fit), 'W_')])]
        colnames(stat) <- sapply(colnames(stat), function(x){str_replace(x,"W_","")})
        
        p.val <- fit[,c("taxon", colnames(fit)[startsWith(colnames(fit), 'p_')])]
        colnames(p.val) <- sapply(colnames(p.val), function(x){str_replace(x,"p_","")})
    
        coeff <- gather(coeff %>%   rename("taxa"="taxon"), "variable", "coefficient", -taxa)
        se <-    gather(se %>%      rename("taxa"="taxon"), "variable", "se", -taxa)
        stat <-  gather(stat %>%    rename("taxa"="taxon"), "variable", "stat", -taxa)
        p.val <- gather(p.val %>%   rename("taxa"="taxon"), "variable", "pvalue", -taxa)
        
        res.full <- merge(merge(merge(coeff, se, by=c("taxa","variable")), stat, by=c("taxa","variable")), p.val, by=c("taxa","variable"))
        names(res.full)[names(res.full)=="coefficient"] <- "effectsize"
        names(res.full)[names(res.full)=="variable"] <- "variable.mda"
        
        res <- mda.common_do(D, f_idx, res.full, "ancombc2")
        res

    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)
    
}
