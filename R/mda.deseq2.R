###############################################################################
# Run DESEQ2


mda.deseq2 <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages(require("DESeq2"))
    suppressPackageStartupMessages(require(tidyr))
    suppressPackageStartupMessages(require(tibble))


    do <- function(f_idx){

        f <- D$formula$norand[[f_idx]]
        mainvar <- D$formula$main_var[f_idx]
        
        if ( (length(D$formula$rand_intercept[[f_idx]]) + length(D$formula$rand_slope[[f_idx]])) > 0 ){
            message(paste0(c("[MDA] mda.deseq2: Formula ", f_idx, " contains random effects. DESeq2 can not handle random effects. Run will continue without random effects.")))
        }
        
        # We need to remove na rows
        terms <- labels(terms(f))
        variables <- terms[!grepl(':', terms)]
        meta_data.nona <- na.omit(D$meta_data[,variables,drop=FALSE])
        count_data.nona <- D$count_data[rownames(meta_data.nona),]


        dds_res <- mda.cache_load_or_run_save(D, "deseq2", f, 
                    {
                    dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(count_data.nona),
                                                          colData = meta_data.nona,
                                                          design = f)
                    DESeq2::DESeq(dds, sfType = "poscounts")
                    })

        
        res.full <- DESeq2::results(dds_res, name=DESeq2::resultsNames(dds_res)[startsWith(DESeq2::resultsNames(dds_res), mainvar)],
                        tidy=T, format="DataFrame")
        names(res.full)[names(res.full)=="row"] <- "taxa"
        res.full$variable <- rep(mainvar, dim(res.full)[1])
        res.full$formula <- rep(mda.deparse(f), dim(res.full)[1])
        res.full$method <- rep("DESeq2", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

        res <- res.full
        res <- res[res$taxa %in% D$nonrare,]
        res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")

        return(list(res=res, res.full=res.full))
    }

    R <- lapply(1:length(D$formula$main_var), do)

    res <- bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")

    names(res)[names(res)=="log2FoldChange"] <- "effectsize"
    names(res)[names(res)=="lfcSE"] <- "se"

    res.full <- bind_rows(lapply(R, function(x){x$res.full}))

    ###############################################################################
    # Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}