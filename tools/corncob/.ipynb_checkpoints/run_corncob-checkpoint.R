###############################################################################
# DESEQ2 run script

library(tidyverse)
library("DESeq2")

source("../common.R")

###############################################################################
# Read input variables

args = commandArgs(trailingOnly=TRUE)

D <- mda.load(args)

###############################################################################
# Run DESEQ2


dds <- DESeq2::DESeqDataSetFromMatrix(countData = D$count_data,
                                      colData = D$meta_data,
                                      design = D$form_norand)
dds_res <- DESeq2::DESeq(dds, sfType = "poscounts")

res <- results(dds_res, tidy=T, format="DataFrame")

rownames(res) <- res$row
res <- res[,-1]

write.table(res, file=args[3], quote=FALSE, sep="\t", col.names = NA)



###############################################################################
# Output

coeff <- gather(as.data.frame(fit$coefficients) %>% rownames_to_column('taxa'), "factor", "coefficient", 2:4)
stdev <- gather(as.data.frame(fit$stdev.unscaled) %>% rownames_to_column('taxa'), "factor", "stdev", 2:4)
p.val <- gather(as.data.frame(fit$p.value) %>% rownames_to_column('taxa'), "factor", "p.value", 2:4)

res <- merge(merge(coeff, stdev, by=c("taxa","factor")), p.val, by=c("taxa","factor"))

write_ctv(res, outfile)