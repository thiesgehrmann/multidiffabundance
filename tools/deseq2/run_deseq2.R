###############################################################################
# DESEQ2 run script

library(tidyverse)
library("DESeq2")

source(paste0(c(dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])), '../common.R'), collapse="/"))


###############################################################################
# Read input variables

args = commandArgs(trailingOnly=TRUE)
#args <- c("~/repos//UA_isala/flow1_redo_revision/data/genus_level_counts.tsv",
#          '~/repos//UA_isala/flow1_redo_revision/data/metadata.tsv',
#          '~/repos//multidiffabundance/tools/list_of_functions.txt',
#          'output.tsv')

D <- mda.load(args)

###############################################################################
# Run DESEQ2

do <- function(f_idx){
    
    f <- D$formula$norand[[f_idx]]
    print(f)
    mainvar <- D$formula$main_var[f_idx]

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = t(D$count_data),
                                          colData = D$meta_data,
                                          design = f)
    dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")

    res.full <- DESeq2::results(dds_res, name=DESeq2::resultsNames(dds_res)[startsWith(DESeq2::resultsNames(dds_res), mainvar)],
                           tidy=T, format="DataFrame")

    names(res.full)[names(res.full)=="row"] <- "taxa"
    res.full$determinant <- rep(mainvar, dim(res.full)[1])
    res.full$formula <- rep(format(f), dim(res.full)[1])
    res.full$method <- rep("DESeq2", dim(res.full)[1])
    res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
    res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

    res <- res.full
    res <- res[res$taxa %in% D$nonrare,]
    res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
    
    return(list(res=res, res.full=res.full))
}

R <- lapply(1:length(D$formula$main_var), do)

res <- bind_rows(lapply(R, function(x){x$res}))
res$qvalue <- p.adjust(res$pvalue, "fdr")

names(res)[names(res)=="log2FoldChange"] <- "effectsize"


res.full <- bind_rows(lapply(R, function(x){x$res.full}))

write_tsv(res, D$outfile)
write_tsv(res.full, paste0(c(D$outfile, ".full.tsv"), collapse=""))



###############################################################################
# Output
