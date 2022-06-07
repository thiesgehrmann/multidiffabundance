###############################################################################
# LIMMA run script

library(tidyverse)
library("edgeR")

source(paste0(c(dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])), '../common.R'), collapse="/"))


###############################################################################
# Read input variables

#args <- c("~/repos//UA_isala/flow1_redo_revision/data/genus_level_counts.tsv",
#          '~/repos//UA_isala/flow1_redo_revision/data/metadata.tsv',
#          '~/repos//multidiffabundance/tools/list_of_functions.txt',
#          'output.tsv')
args = commandArgs(trailingOnly=TRUE)

D <- mda.load(args)

###############################################################################
#


DGE_LIST <- DGEList(t(D$count_data))

### check if upper quartile method works for selecting reference
Upper_Quartile_norm_test <- calcNormFactors(DGE_LIST, method="upperquartile")

summary_upper_quartile <- summary(Upper_Quartile_norm_test$samples$norm.factors)[3]
if(is.na(summary_upper_quartile) | is.infinite(summary_upper_quartile)){
    message("Upper Quartile reference selection failed will use find sample with largest sqrt(read_depth) to use as reference")
    Ref_col <- which.max(colSums(sqrt(t(D$count_data))))
    DGE_LIST_Norm <- calcNormFactors(DGE_LIST, method = "TMM", refColumn = Ref_col)
    fileConn<-file(args[[4]])
    writeLines(c("Used max square root read depth to determine reference sample"), fileConn)
    close(fileConn)

}else{
    DGE_LIST_Norm <- calcNormFactors(DGE_LIST, method="TMM")
}

do <- function(f_idx){
    
    f <- D$formula$norand[[f_idx]]
    print(f)
    mainvar <- D$formula$main_var[f_idx]

    ## make matrix for testing
    mm <- model.matrix(f, D$meta_data)

    voomfit <- voom(DGE_LIST_Norm, mm, plot=F)

    fit <- lmFit(voomfit, mm)
    fit <- eBayes(fit)

    # Gather output

    coeff <- gather(as.data.frame(fit$coefficients) %>% rownames_to_column('taxa'), "determinant", "coefficient", 2:(dim(as.data.frame(fit$coefficients))[2]+1))
    stdev <- gather(as.data.frame(fit$stdev.unscaled) %>% rownames_to_column('taxa'), "determinant", "stdev", 2:(dim(as.data.frame(fit$stdev.unscaled))[2]+1))
    p.val <- gather(as.data.frame(fit$p.value) %>% rownames_to_column('taxa'), "determinant", "pvalue", 2:(dim(as.data.frame(fit$p.value))[2]+1))

    res.full <- merge(merge(coeff, stdev, by=c("taxa","determinant")), p.val, by=c("taxa","determinant"))
    res.full$formula <- rep(format(f), dim(res.full)[1])
    res.full$method <- rep("limmaVoom", dim(res.full)[1])
    res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
    res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])
    
    # Select only the relevant determinant & taxa
    res <- res.full
    res <- res[startsWith(res$determinant, mainvar),]
    res <- res[res$taxa %in% D$nonrare,]

    res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
    res$determinant <- rep(mainvar, dim(res)[1])
    names(res)[names(res)=="coefficient"] <- "effectsize"

    return(list(res=res, res.full=res.full))
}

R <- lapply(1:length(D$formula$main_var), do)

###############################################################################
# Output

res <- bind_rows(lapply(R, function(x){x$res}))
res$qvalue <- p.adjust(res$pvalue, "fdr")
res.full <- bind_rows(lapply(R, function(x){x$res.full}))

write_tsv(res, D$outfile)
write_tsv(res.full, paste0(c(D$outfile, ".full.tsv"), collapse=""))