mda.ancombc2 <- function(mda.D, ...){
    suppressPackageStartupMessages({
        require("ANCOMBC")
        require(tidyr)
        require(dplyr)
        require(tibble)
        require(stringr)})
    
    D <- mda.D

    OTU <- phyloseq::otu_table(t(D$count_data), taxa_are_rows = T)
    sampledata <- phyloseq::sample_data(D$meta_data, errorIfNULL = F)
    phylo <- phyloseq::merge_phyloseq(OTU, sampledata)

    do <- function(f_idx){
        f <- D$formula$norand[[f_idx]]
        
        f.ancombc.fix <- strsplit(mda.deparse(f),"~")[[1]][2]
        f.ancombc.rand <- paste0(c(D$formula$rand_intercept[[f_idx]], D$formula$rand_slope[[f_idx]]), collapse="+")
        f.ancombc.rand <- if (f.ancombc.rand == ''){NULL} else {f.ancombc.rand}

        mainvar <- D$formula$main_var[f_idx]

        out <- mda.cache_load_or_run_save(D, "ancombc", f, 
                   ANCOMBC::ancombc2(data = phylo, fix_formula = f.ancombc.fix, rand_formula = f.ancombc.rand, 
                           p_adj_method = "holm", prv_cut=0, lib_cut = 1000, 
                           struc_zero = FALSE, global = FALSE, alpha = 0.05) )

        fit = out$res
        
        coeff <- fit[,c("taxon", colnames(fit)[startsWith(colnames(fit), 'lfc_')])]
        colnames(coeff) <- sapply(colnames(coeff), function(x){str_replace(x,"lfc_","")})
        
        se <- fit[,c("taxon", colnames(fit)[startsWith(colnames(fit), 'se_')])]
        colnames(se) <- sapply(colnames(se), function(x){str_replace(x,"se_","")})
        
        stat <- fit[,c("taxon", colnames(fit)[startsWith(colnames(fit), 'W_')])]
        colnames(stat) <- sapply(colnames(stat), function(x){str_replace(x,"W_","")})
        
        p.val <- fit[,c("taxon", colnames(fit)[startsWith(colnames(fit), 'p_')])]
        colnames(p.val) <- sapply(colnames(p.val), function(x){str_replace(x,"p_","")})

        coeff <- gather(coeff %>%   rename("taxa"="taxon"), "variable", "coefficient", -taxa)
        se <-    gather(se %>%    rename("taxa"="taxon"), "variable", "se", -taxa)
        stat <-  gather(stat %>%     rename("taxa"="taxon"), "variable", "stat", -taxa)
        p.val <- gather(p.val %>% rename("taxa"="taxon"), "variable", "pvalue", -taxa)
        #colnames(coeff) <- colnames(coeff) <- colnames(coeff) <- colnames(coeff) <- c("taxa","variable","coefficient")
        
        #coeff <- gather(as.data.frame(fit$lfc) %>% rownames_to_column('taxa'), "variable", "coefficient", 2:(dim(as.data.frame(fit$lfc))[2]+1))
        #se <- gather(as.data.frame(fit$se) %>% rownames_to_column('taxa'), "variable", "se", 2:(dim(as.data.frame(fit$se))[2]+1))
        #stat <- gather(as.data.frame(fit$W) %>% rownames_to_column('taxa'), "variable", "stat", 2:(dim(as.data.frame(fit$W))[2]+1))
        #p.val <- gather(as.data.frame(fit$p_val) %>% rownames_to_column('taxa'), "variable", "pvalue", 2:(dim(as.data.frame(fit$p_val))[2]+1))

        #return(coeff)
        res.full <- merge(merge(merge(coeff, se, by=c("taxa","variable")), stat, by=c("taxa","variable")), p.val, by=c("taxa","variable"))
        res.full$formula <- rep(mda.deparse(f), dim(res.full)[1])
        res.full$method <- rep("ancombc", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

        res <- res.full
        res <- res[startsWith(res$variable, mainvar),]
        res <- res[res$taxa %in% D$nonrare,]

        res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
        res$variable <- rep(mainvar, dim(res)[1])
        names(res)[names(res)=="coefficient"] <- "effectsize"

        return(list(res=res, res.full=res.full))
    }

    R <- lapply(1:length(D$formula$main_var), do)
    #return(R)

    res <- bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    res.full <- bind_rows(lapply(R, function(x){x$res.full}))

    ###############################################################################
    # Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}