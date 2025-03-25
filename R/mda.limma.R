#' @export
mda.limma <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages({
        require("edgeR")
        require(tidyr)
        require(tibble)
        require(dplyr)})

    DGE_LIST_Norm <- mda.cache_load_or_run_save(D, NULL, "limma_dgelistnorm", {
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
        DGE_LIST_Norm
    })
    

    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn
        
        if ( length(fdata$parts.random.slope) > 0 ){
            message(paste0(c("[MDA] mda.limma: Formula on ", f_idx, " contains random slope effects. limma can not handle random slopes.")))
            return(mda.common_do(D, f_idx, mda.empty_output(D, f_idx, "Formula incompatible with limma analysis (random slope specified)"), "limma", skip_taxa_sel=TRUE))
        }
        
        if ( length(fdata$parts.random.intercept) > 1 ){
            message(paste0(c("[MDA] mda.limma: Formula on ", f_idx, " contains more than one random intercept effect. limma can only handle one random intercept")))
            return(mda.common_do(D, f_idx, mda.empty_output(D, f_idx, "Formula incompatible with limma analysis (>1 random intercept specified)"), "limma", skip_taxa_sel=TRUE))
        }
        
        block <- NULL
        if ( length(fdata$parts.random.intercept) == 1 ){
            rblock <- fdata$parts.random.intercept[1]
            block <- trimws(unlist(strsplit(rblock, split="\\|"))[[2]])
            block <- fdata$data[,block]
        }

        fit <- mda.cache_load_or_run_save(D, f_idx, "limma", 
                   {
                    mm <- fdata$data[,fdata$parts.fixed,drop=FALSE]
                    mm <- mm[complete.cases(mm), ,drop=FALSE] # Remove rows with NA... (Note: Do this in the formulas prep, or here? all tools should drop the NA rows anyway, right? Look into this...)

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
        names(res.full)[names(res.full)=="coefficient"] <- "effectsize"
        names(res.full)[names(res.full)=="stdev"] <- "se"
        names(res.full)[names(res.full)=="variable"] <- "variable.mda"
        
        res <- mda.common_do(D, f_idx, res.full, "limma", skip_taxa_sel=FALSE)

        res$res.full$se <- res$res.full$se / sqrt(as.numeric(res$res.full$n))
        res$res$se <- res$res$se / sqrt(as.numeric(res$res$n))
        res

    }

    R <- lapply(1:length(D$formula), do)

    mda.common_output(R)
}
