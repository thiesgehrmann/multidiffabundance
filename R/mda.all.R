
mda.all <- function(mda.D, alpha=FALSE, beta=FALSE, group=FALSE, continuous=FALSE,
                    aldex2=TRUE, ancombc2=TRUE, corncob=FALSE, deseq2=TRUE,
                    limma=TRUE, lmclr=TRUE, maaslin2=TRUE, zicoseq=TRUE, ...){
    suppressPackageStartupMessages(require(dplyr))
    D <- mda.D
    
    functions <- c(mda.alpha, mda.beta, mda.group, mda.continuous, mda.aldex2, mda.ancombc2, mda.corncob, mda.deseq2, mda.limma, mda.lmclr, mda.maaslin2, mda.zicoseq)
    selected <-  c(    alpha,     beta,     group,     continuous,     aldex2,     ancombc2,      corncob,     deseq2,     limma,     lmclr,     maaslin2,     zicoseq)
    
    R <- lapply(functions[selected], function(x){ x(D, ...) })
    res <- bind_rows(lapply(R, function(x){x$res}))
    res.full <- bind_rows(lapply(R, function(x){x$res.full}))
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))

}