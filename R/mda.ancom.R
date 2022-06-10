###############################################################################
# Run ANCOMBC

mda.ancombc <- function(mda.D){
    D <- mda.D
    require("ANCOMBC")

    OTU <- phyloseq::otu_table(t(D$count_data), taxa_are_rows = T)
    sampledata <- phyloseq::sample_data(D$meta_data, errorIfNULL = F)
    phylo <- phyloseq::merge_phyloseq(OTU, sampledata)

    do <- function(f_idx){
        f <- strsplit(format(D$formula$norand[[f_idx]]),"~")[[1]][2]
        print(f)
        mainvar <- D$formula$main_var[f_idx]

        out <- mda.cache_load_or_run_save(D$cacheprefix, "ancombc", f, 
                   ancombc(phyloseq = phylo, formula = f, 
                           p_adj_method = "holm", zero_cut = 0.90, lib_cut = 1000, 
                           struc_zero = FALSE, global = FALSE, neg_lb = TRUE, tol = 1e-5, 
                           max_iter = 100, conserve = TRUE, alpha = 0.05) )

        fit = out$res

        coeff <- gather(as.data.frame(fit$beta) %>% rownames_to_column('taxa'), "variable", "coefficient", 2:(dim(as.data.frame(fit$beta))[2]+1))
        se <- gather(as.data.frame(fit$se) %>% rownames_to_column('taxa'), "variable", "se", 2:(dim(as.data.frame(fit$se))[2]+1))
        stat <- gather(as.data.frame(fit$W) %>% rownames_to_column('taxa'), "variable", "stat", 2:(dim(as.data.frame(fit$W))[2]+1))
        p.val <- gather(as.data.frame(fit$p_val) %>% rownames_to_column('taxa'), "variable", "pvalue", 2:(dim(as.data.frame(fit$p_val))[2]+1))

        res.full <- merge(merge(merge(coeff, se, by=c("taxa","variable")), stat, by=c("taxa","variable")), p.val, by=c("taxa","variable"))
        res.full$formula <- rep(format(f), dim(res.full)[1])
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

    res <- bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    res.full <- bind_rows(lapply(R, function(x){x$res.full}))

    ###############################################################################
    # Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    
    return(list(res=res, res.full=res.full))
}