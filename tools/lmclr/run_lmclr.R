###############################################################################
# LIMMA run script

library(tidyverse)

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
# LMCLR

lmclr <- function(count_data, meta_data, formula, mainvar, taxa=NULL){
    taxa <- if (is.null(taxa)) colnames(count_data) else taxa
    clr_data <- mda.clr(mda.pseudocount(count_data))
    
    f <- update(formula, clrtaxa ~ .)
    
    res <- lapply(taxa, function(t){
        meta_data$clrtaxa <- clr_data[,t]
        fit <- lm(f, data=meta_data, na.action = 'na.exclude')
        s <- as.data.frame(coefficients(summary(fit)))
        s$taxa <- rep(t, dim(s)[1])
        s <- s %>% rownames_to_column("determinant")
        s
    })
    res <- bind_rows(res)

    names(res)[names(res)=="coefficient"] <- "effectsize"
    names(res)[names(res)=="Std. Error"] <- "se"
    names(res)[names(res)=="t value"] <- "stat"
    names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"
    
    res
}

do <- function(f_idx){
    
    f <- D$formula$norand[[f_idx]]
    print(f)
    mainvar <- D$formula$main_var[f_idx]

    res.full <- lmclr(D$count_data, D$meta_data, f, mainvar, D$nonrare)
    res.full$formula <- rep(format(f), dim(res.full)[1])
    res.full$method <- rep("lmclr", dim(res.full)[1])
    res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
    res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])
    
    # Select only the relevant determinant ( taxa are selected already in lmclr )
    res <- res.full
    res <- res[startsWith(res$determinant, mainvar),]


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