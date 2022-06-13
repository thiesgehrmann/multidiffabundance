
mda.all <- function(mda.D, aldex2=TRUE, ancombc=TRUE, corncob=TRUE, deseq2=TRUE, limma=TRUE, lmclr=TRUE, maaslin2=TRUE){
    require(dplyr)
    D <- mda.D
    
    functions <- c(mda.aldex2, mda.ancombc, mda.corncob, mda.deseq2, mda.limma, mda.lmclr, mda.maaslin2)
    selected <- c(aldex2, ancombc, corncob, deseq2, limma, lmclr, maaslin2)
    
    R <- lapply(functions[selected], function(x){ x(D) })
    res <- bind_rows(lapply(R, function(x){x$res}))
    res.full <- lapply(R, function(x){x$res.full})
    
    return(list(res=res, res.full=res.full))

}