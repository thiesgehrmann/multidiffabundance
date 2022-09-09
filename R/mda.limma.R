###############################################################################
# RUN LIMMA

mda.limma <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages(require("edgeR"))
    suppressPackageStartupMessages(require(tidyr))
    suppressPackageStartupMessages(require(tibble))
    suppressPackageStartupMessages(require(dplyr))

    DGE_LIST <- DGEList(t(D$count_data))

    ### check if upper quartile method works for selecting reference
    Upper_Quartile_norm_test <- calcNormFactors(DGE_LIST, method="upperquartile")

    summary_upper_quartile <- summary(Upper_Quartile_norm_test$samples$norm.factors)[3]
    if(is.na(summary_upper_quartile) | is.infinite(summary_upper_quartile)){
        message("Upper Quartile reference selection failed will use find sample with largest sqrt(read_depth) to use as reference")
        Ref_col <- which.max(colSums(sqrt(t(D$count_data))))
        DGE_LIST_Norm <- calcNormFactors(DGE_LIST, method = "TMM", refColumn = Ref_col)
        message("Used max square root read depth to determine reference sample")

    }else{
        DGE_LIST_Norm <- calcNormFactors(DGE_LIST, method="TMM")
    }

    do <- function(f_idx){

        f <- D$formula$norand[[f_idx]]
        f.cache <- f

        mainvar <- D$formula$main_var[f_idx]
        
        if ( length(D$formula$rand_slope[[f_idx]]) > 0 ){
            message(paste0(c("[MDA] mda.limma: Formula on ", f_idx, " contains random slope effects. limma can not handle random slopes. Run will continue without random slopes")))
        }
        
        block <- NULL
        if ( length(D$formula$rand_intercept[[f_idx]]) > 0 ){
            rblock <- D$formula$rand_intercept[[f_idx]][1]
            block <- trimws(gsub(")", "", unlist(strsplit(rblock, split="\\|"))[[2]]))
            f.cache <- update.formula(f.cache, paste0(c("~.+", paste0(rblock, collapse="+")), collapse=""))
            if (length(D$formula$rand_intercept[[f_idx]]) > 1){
                message(paste0(c("[MDA] mda.limma: Formula on ", f_idx, " contains more than one random intercept effect. limma can only handle one random effect. Will continue with (1|", block, ").")))
            }
            block <- D$meta_data[,block]
        }
        

        ## make matrix for testing
        fit <- mda.cache_load_or_run_save(D, "limma", f.cache, 
                   {
                    mm <- model.matrix(f, D$meta_data)
                    subset_dgelist <- DGE_LIST_Norm[, rownames(DGE_LIST_Norm$samples) %in% rownames(mm)]
                    voomfit <- voom(subset_dgelist, mm, plot=FALSE)
                    fit <- if (is.null(block)) {
                        lmFit(voomfit, mm)
                    } else {
                        # With explanations on 
                        # https://bioconductor.riken.jp/packages/3.9/bioc/vignettes/variancePartition/inst/doc/dream.html
                        subset_block <- block[rownames(DGE_LIST_Norm$samples) %in% rownames(mm)]
                        dupcor <- duplicateCorrelation(voomfit,mm,block=subset_block)
                        vobj = voom( subset_dgelist, mm, plot=FALSE, block=subset_block, correlation=dupcor$consensus)
                        dupcor <- duplicateCorrelation(vobj, mm, block=subset_block)
                        lmFit(vobj, mm, block=subset_block, correlation=dupcor$consensus)
                    }
                    eBayes( fit )
                   } )

        # Gather output

        coeff <- gather(as.data.frame(fit$coefficients) %>% rownames_to_column('taxa'), "variable", "coefficient", 2:(dim(as.data.frame(fit$coefficients))[2]+1))
        stdev <- gather(as.data.frame(fit$stdev.unscaled) %>% rownames_to_column('taxa'), "variable", "stdev", 2:(dim(as.data.frame(fit$stdev.unscaled))[2]+1))
        p.val <- gather(as.data.frame(fit$p.value) %>% rownames_to_column('taxa'), "variable", "pvalue", 2:(dim(as.data.frame(fit$p.value))[2]+1))
        stat <- gather(as.data.frame(fit$t) %>% rownames_to_column('taxa'), "variable", "stat", 2:(dim(as.data.frame(fit$t))[2]+1))


        res.full <- merge(merge(merge(coeff, stdev, by=c("taxa","variable")), p.val, by=c("taxa","variable")), stat, by=c("taxa","variable"))
        res.full$formula <- rep(mda.deparse(f.cache), dim(res.full)[1])
        res.full$method <- rep("limma", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])
        res.full$se <- res.full$stdev / sqrt(res.full$n)

        # Select only the relevant variable & taxa
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
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}