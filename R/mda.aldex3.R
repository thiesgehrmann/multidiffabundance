############################################
# ALDEx3 function
mda.aldex3 <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages({
        require(ALDEx3)
        require(dplyr)
        require(tibble)})

    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]

        abun <- as.data.frame(D$count_data)
        meta <- fdata$data
        form <- fdata$fn

        nonrare <- abun[,D$nonrare]
        other <- rowSums(abun) - rowSums(nonrare)
        nonrare$other <- other

        data.in <- as.data.frame(t(nonrare))

        # OK so this is some strange behaviour of aldex3, apparently it can't handle NAs in the metadata
        # So cleaning this up:
        meta.nona <- tidyr::drop_na(meta)
        data.in.nona <- data.in[,rownames(meta.nona)]

        r <- mda.trycatchempty(D, f_idx, {
            mda.cache_load_or_run_save(D, f_idx, "aldex3", {
                aldex.fit <- aldex(data.in.nona,
                                   form,
                                   meta.nona,
                                   nsample=2000,
                                   scale=clr.sm,   # CLR assumption
                                   gamma=0,        # Gamma=0 no scale uncertainty
                                   return.pars=c("X", "estimate", "std.error", "p.val",
                                         "p.val.adj", "logComp", "logScale"))
                })
            
            }, taxa=D$nonrare)
            
        res <- if (r$error){
            mda.message(r$message, type="error")
            r$response
        } else {
            aldex.fit <- r$response
            
            effectsize <- reshape2::melt(aldex.fit$estimate, value.name="estimate")
            se         <- reshape2::melt(aldex.fit$std.error, value.name="se")
            pvalue     <- reshape2::melt(aldex.fit$p.val, value.name="pvalue")

            if (all(effectsize$Var2 == pvalue$Var2)) {
              merged <- effectsize
              merged$se <- se$se
              merged$pvalue <- pvalue$pvalue
                
            } else {
              merged <- merge(effectsize, merge(se, pvalue, by=c("Var1",'Var2')), by=c("Var1",'Var2'))
            }
            
            
            colnames(merged) <- c("variable.mda", 'taxa', "effectsize", 'se', 'pvalue')
            merged
        }
            
        res.full <- mda.common_do(D, f_idx, as.data.frame(res), "aldex3", skip_taxa_sel=FALSE)
        res.full
    }

    R <- lapply(1:length(D$formula), do)
    mda.common_output(R)
}
