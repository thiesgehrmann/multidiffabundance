###############################################################################
# ALDEX2 run script

library(readr)
library(tidyverse)
library("ALDEx2")

source(paste0(c(dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])), '../common.R'), collapse="/"))

###############################################################################
# Read input variables
args <- c("~/repos//UA_isala/flow1_redo_revision/data/genus_level_counts.tsv",
          '~/repos//UA_isala/flow1_redo_revision/data/metadata.tsv',
          '~/repos//multidiffabundance/tools/list_of_functions.txt',
          'output.tsv')
args = commandArgs(trailingOnly=TRUE)

D <- mda.load(args)

###############################################################################
# ALDEx2

do <- function(f_idx){
    f <- D$formula$norand[[f_idx]]
    print(f)
    mainvar <- D$formula$main_var[f_idx]
    
    
    out <- mda.cache_load_or_run_save(D$cacheprefix, "aldex2", f, {
        model_matrix <- model.matrix(f, D$meta_data)
        aldex(t(D$count_data), model_matrix, denom = "all", test="glm")
    })
    
    res.full <- gather(as.data.frame(out) %>% rownames_to_column('taxa'), "measure", "value", 2:(dim(as.data.frame(out))[2]+1))
    
    clean.feature <- function(v){
        ret <- if (endsWith(v, ".Estimate")) {
                    "effectsize"
                } else if (endsWith(v, ".Std..Error")) {
                    "se"
                } else if (endsWith(v, ".t.value")) {
                    "stat"
                } else if (endsWith(v, ".Pr...t..")) {
                    "pvalue"
                } else if (endsWith(v, ".Pr...t...BH")) {
                    "qvalue"
                }
        ret
    }
    res.full$feature <- unlist(lapply(res.full$measure, clean.feature))
    
    clean.determinant <- function(v){
        gsub("^model.", "", gsub(".Estimate", "", gsub(".Pr...t..", "", gsub(".Pr...t...BH", "", gsub(".t.value", "", gsub(".Std..Error", "", v))))))
    }
    res.full$determinant <- unlist(lapply(res.full$measure, clean.determinant))
    
    res.full <- pivot_wider(res.full, id_cols=c("taxa", "determinant"), names_from=feature, values_from=value)
    res.full$formula <- rep(format(f), dim(res.full)[1])
    res.full$method <- rep("aldex2", dim(res.full)[1])
    res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
    res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])
    
    
    res <- res.full[startsWith(res.full$determinant, mainvar),]
    res$determinant <- rep(mainvar, dim(res)[1])
    res <- res[res$taxa %in% D$nonrare,]

    res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")

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