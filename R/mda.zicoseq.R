###############################################################################
#' ZICOSEQ
#' @export
mda.zicoseq <- function(mda.D, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(GUniFrac)
        require(reshape2)
        require(tibble)
        require(dplyr)})

    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn
        
        if ( length(fdata$parts.random) > 0 ){
            message(paste0(c("[MDA] mda.zicoseq: Formula on ", f_idx, " contains random effects. ZicoSeq can not handle random effects")))
            return(mda.common_do(D, f_idx, mda.empty_output(D, f_idx, "Formula incompatible with ZicoSeq analysis (random effect specified)", taxa=D$nonrare), "zicoseq", skip_taxa_sel=FALSE))
        }
        
        sel_taxa <- colnames(D$count_data)[apply(D$count_data, 2, sd) != 0]
        
        response <- mda.trycatchempty(D, f_idx, {
            mda.cache_load_or_run_save(D, f_idx, "zicoseq", {
                mda.empty_output(D, 1, taxa=D$nonrare, comment="zicoseq has a stupid behaviour")

                if (length(fdata$parts.fixed) == 1) {
                    adjust <- NULL
                } else {
                    adjust <- fdata$parts.fixed[-1]
                }
                results <- rlang::exec(ZicoSeq,
                    meta.dat = fdata$data,
                    feature.dat = t(D$count_data[, sel_taxa]),
                    grp.name = fdata$parts.fixed[1],
                    adj.name = adjust,
                    feature.dat.type = "count",
                    # Filter to remove rare taxa
                    min.prop = D$nonrare.pct,
                    # Winsorization to replace outliers
                    # Posterior sampling
                    # Use the square-root transformation
                    # Permutation-based multiple testing correction
                    # Reference-based multiple stage normalization
                    # Family-wise error rate control
                    !!!D$args$zicoseq,
                    ...
                    )

                res <- melt(as.data.frame(t(results$coef.list[[1]])) %>% rownames_to_column('taxa'), variable.name="variable.mda", value.name="effectsize", id.vars="taxa")

                res.pval <- melt(results$p.raw) %>% rownames_to_column('taxa')
                colnames(res.pval) <- c("taxa", 'pvalue')
                res.R2 <- melt(results$R2)[,c('Var1','value')]
                colnames(res.R2) <- c("taxa", 'stat')
                res.other <- merge(res.pval, res.R2, by="taxa")
                res.other$variable.mda <- results$grp.name

                results <- merge(res, res.other, by=c('taxa','variable.mda'))
                results
                }, order_invariant=FALSE)
            }, taxa=D$nonrare)
        v <- mda.common_do(D, f_idx, response$response, "zicoseq", skip_taxa_sel=FALSE)
        v
    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)
}
