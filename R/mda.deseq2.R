###############################################################################
#' Run DESEQ2
#' @export
mda.deseq2 <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages({
        require("DESeq2")
        require(tidyr)
        require(tibble)})

    do <- function(f_idx){

        fdata <- D$formula[[f_idx]]
        f <- fdata$fn
        
        if ( length(fdata$parts.random) > 0 ){
            message(paste0(c("[MDA] mda.deseq2: Formula on ", f_idx, " contains random effects. DESeq2 can not handle random effects")))
            return(mda.common_do(D, f_idx, mda.empty_output(D, f_idx, "Formula incompatible with DESeq2 analysis (random effect specified)"), "deseq2", skip_taxa_sel=TRUE))
        }
        
        # We need to remove na rows

        meta_data.nona <- na.omit(fdata$data)
        count_data.nona <- D$count_data[rownames(meta_data.nona),]

        r <- mda.trycatchempty(D, f_idx, {
                mda.cache_load_or_run_save(D, f_idx, "deseq2", 
                        {
                        dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(count_data.nona),
                                                              colData = meta_data.nona,
                                                              design = f)
                        rlang::exec(DESeq2::DESeq, dds, !!!D$args$deseq2)
                        })
            }, taxa=D$nonrare)
            
        if (r$error){
            return(mda.common_do(D, f_idx, r$response, "deseq2", skip_taxa_sel = FALSE))
        }

        dds_res = r$response


        res.full <- dplyr::bind_rows(lapply(DESeq2::resultsNames(dds_res), function(name){
            v <- DESeq2::results(dds_res, name=name, tidy=T, format="DataFrame")
            v$variable.mda <- name
            v}))
        names(res.full)[names(res.full)=="row"] <- "taxa"
        names(res.full)[names(res.full)=="log2FoldChange"] <- "effectsize"
        names(res.full)[names(res.full)=="lfcSE"] <- "se"

        mda.common_do(D, f_idx, res.full, "deseq2", skip_taxa_sel=FALSE)
    }

    R <- lapply(1:length(D$formula), do)

    mda.common_output(R)
}
