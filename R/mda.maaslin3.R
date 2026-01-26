###############################################################################
# Maaslin3

mda.maaslin3 <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages({
        require(maaslin3)
        require(dplyr)
        require(tibble)})

    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]

        ### NEED TO VERIFY IF THIS WORKS. FOR NOW TEST IS EXCLUDED.
        #if ( length(fdata$parts.random.slope) > 0 ){
        #    message(paste0(c("[MDA] mda.maaslin3: Formula ", f_idx, " contains random slope effects. Maaslin2 can not handle random slopes.")))
        #    return(mda.common_do(D, f_idx, mda.empty_output(D, f_idx, "Formula incompatible with maaslin2 analysis (random slope specified)"), "maaslin3", skip_taxa_sel=TRUE))
        #}

        r.ab <- D$count_data / rowSums(D$count_data)
        meta <- fdata$data
        form <- fdata$fn
        out.dir <- paste0(c(D$outprefix, "/maaslin3.output.folder"), collapse="")
        
        fit_out <- mda.cache_load_or_run_save(D, f_idx, "maaslin3",
                    maaslin3(
                        input_data = r.ab,
                        input_metadata = meta,
                        output = out.dir,
                        formula = form,
                        normalization = 'TSS',
                        transform = 'LOG',
                        augment = TRUE,
                        standardize = TRUE,
                        max_significance = 0.1,
                        median_comparison_abundance = TRUE,
                        median_comparison_prevalence = FALSE,
                        max_pngs = 100,
                        cores = 1,
                        save_models = TRUE))

        res.ab <- fit_out$fit_data_abundance$results[,c('feature','name','coef','stderr','pval_individual','error')]
        res.ab$qvalue <- p.adjust(res.ab$pval_individual, 'fdr')
        colnames(res.ab) <- c("taxa", 'variable.mda', 'effectsize', 'se', 'pvalue', 'comment', 'qvalue')
        
        res.pr <- fit_out$fit_data_prevalence$results[,c('feature','name','coef','stderr','pval_individual','error')]
        res.pr$qvalue <- p.adjust(res.pr$pval_individual, 'fdr')
        colnames(res.pr) <- c("taxa", 'variable.mda', 'effectsize', 'se', 'pvalue', 'comment', 'qvalue')
        
        res.j <- fit_out$fit_data_prevalence$results[,c('feature','name','coef','stderr','pval_joint','error')]
        res.j$qvalue <- p.adjust(res.j$pval_joint, 'fdr')
        colnames(res.j) <- c("taxa", 'variable.mda', 'effectsize', 'se', 'pvalue', 'comment', 'qvalue')
        
        res.full <- rbind(res.ab, res.pr, res.j)
        
        res.full <- as.data.frame(res.full)
        
        res.full.ab <- mda.common_do(D, f_idx, as.data.frame(res.ab), "maaslin3.abundance", skip_taxa_sel=FALSE)
        res.full.pr <- mda.common_do(D, f_idx, as.data.frame(res.pr), "maaslin3.prevalence", skip_taxa_sel=FALSE)
        res.full.j  <- mda.common_do(D, f_idx, as.data.frame(res.j), "maaslin3.joint", skip_taxa_sel=FALSE)
        list(abundance=res.full.ab, prevalence=res.full.pr, joint=res.full.j)
    }

    R <- lapply(1:length(D$formula), do)
    R.ab <- lapply(R, function(x){x[["abundance"]]})
    R.pr <- lapply(R, function(x){x[["prevalence"]]})
    R.j  <- lapply(R, function(x){x[["joint"]]})
    

    co.ab <- mda.common_output(R.ab)
    co.pr <- mda.common_output(R.pr)
    co.j  <- mda.common_output(R.j)

    co <- list(co.ab, co.pr, co.j)
    mda.merge_results(co)
}
