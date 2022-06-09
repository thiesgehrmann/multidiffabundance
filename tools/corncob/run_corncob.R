###############################################################################
# DESEQ2 run script

library(corncob)
library(phyloseq)

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
# Run CORNCOB

OTU <- phyloseq::otu_table(t(D$count_data), taxa_are_rows = T)
sampledata <- phyloseq::sample_data(D$meta_data, errorIfNULL = F)
phylo <- phyloseq::merge_phyloseq(OTU, sampledata)

do <- function(f_idx){
    f <- D$formula$norand[[f_idx]]
    print(f)
    mainvar <- D$formula$main_var[f_idx]

    results <- mda.cache_load_or_run_save(D$cacheprefix, "corncob", f,
                   corncob::differentialTest(formula= f,
                                             phi.formula = f,
                                             phi.formula_null = f,
                                             formula_null = ~ 1,
                                             test="Wald", data=phylo,
                                             boot=F,
                                             fdr_cutoff = 0.05) )

    ram <- results$all_models
    ram <- lapply(1:length(ram), function(t_idx){
        ms <- ram[[t_idx]]
        if (!('coefficients' %in% names(ms))){
            return(NULL)
        }
        ms <- as.data.frame(ms$coefficients)
        ms <- ms[startsWith(rownames(ms), "phi."),]
        ms$taxa <- rep(colnames(D$count_data)[t_idx], dim(ms)[1])
        ms <- ms %>% rownames_to_column('determinant')
        names(ms)[names(ms)=="Estimate"] <- "effectsize"
        names(ms)[names(ms)=="Std. Error"] <- "se"
        names(ms)[names(ms)=="t value"] <- "stat"
        names(ms)[names(ms)=="Pr(>|t|)"] <- "pvalue"
        ms
    })
    
    res.full <- bind_rows(ram)
    res.full$formula <- rep(format(f), dim(res.full)[1])
    res.full$method <- rep("corncob", dim(res.full)[1])
    res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
    res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])
    
    res <- res.full
    res <- res[startsWith(res$determinant, mainvar),]
    res <- res[res$taxa %in% D$nonrare,]

    res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
    res$determinant <- rep(mainvar, dim(res)[1])

    return(list(res=res, res.full=res.full))
    
}

R <- lapply(1:length(D$formula$main_var), do)

res <- bind_rows(lapply(R, function(x){x$res}))
res$qvalue <- p.adjust(res$pvalue, "fdr")
res.full <- bind_rows(lapply(R, function(x){x$res.full}))

###############################################################################
# Output

write_tsv(res, paste0(c(D$outprefix, "results.tsv"), collapse=""))
write_tsv(res.full, paste0(c(D$outprefix, "results.full.tsv"), collapse=""))


