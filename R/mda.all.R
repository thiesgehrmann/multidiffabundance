
mda.all <- function(mda.D, alpha=FALSE, beta=FALSE,
                    aldex2=TRUE, ancombc=TRUE, corncob=FALSE, deseq2=TRUE,
                    limma=TRUE, lmclr=TRUE, maaslin2=TRUE, ...){
    suppressPackageStartupMessages(require(dplyr))
    D <- mda.D
    
    functions <- c(mda.alpha, mda.beta, mda.aldex2, mda.ancombc2, mda.corncob, mda.deseq2, mda.limma, mda.lmclr, mda.maaslin2)
    selected <- c(alpha, beta, aldex2, ancombc, corncob, deseq2, limma, lmclr, maaslin2)
    
    R <- lapply(functions[selected], function(x){ x(D, ...) })
    res <- bind_rows(lapply(R, function(x){x$res}))
    res.full <- lapply(R, function(x){x$res.full})
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))

}