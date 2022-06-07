###############################################################################
# DESEQ2 run script

library(corncob)
library(phyloseq)

source("../common.R")

###############################################################################
# Read input variables

args = commandArgs(trailingOnly=TRUE)

D <- mda.load(args)

###############################################################################
# Run CORNCOB

OTU <- phyloseq::otu_table(D$count_data, taxa_are_rows = T)
sampledata <- phyloseq::sample_data(D$meta_data, errorIfNULL = F)
phylo <- phyloseq::merge_phyloseq(OTU, sampledata)


results <- corncob::differentialTest(formula= D$form_nonrand,
                                     phi.formula = D$form_nonrand,
                                     phi.formula_null = D$form_nonrand,
                                     formula_null = ~ 1,
                                     test="Wald", data=phylo,
                                     boot=F,
                                     fdr_cutoff = 0.05)


###############################################################################
# Output

write.table(results$p_fdr, file=args[[3]], sep="\t", col.names = NA, quote=F)
write.table(results$p, file=paste0(args[[3]], "_uncor", sep=""), sep="\t", col.names = NA, quote=F)



