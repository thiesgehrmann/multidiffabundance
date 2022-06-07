###############################################################################
# ALDEX2 run script

library(readr)
library(tidyverse)
library("ALDEx2")

source(paste0(c(dirname(sub("--file=", "", commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))])), '../common.R'), collapse="/"))

###############################################################################
# Read input variables

args = commandArgs(trailingOnly=TRUE)
args <- c("~/repos//UA_isala/flow1_redo_revision/data/genus_level_counts.tsv",
          '~/repos//UA_isala/flow1_redo_revision/data/metadata.tsv',
          '~/repos//multidiffabundance/tools/list_of_functions.txt',
          'output.tsv')

D <- mda.load(args)

###############################################################################
# ALDEx2

do <- function(f_idx){
    f <- D$formula$norand[[f_idx]]
    print(f)
    mainvar <- D$formula$main_var[f_idx]
    
    model_matrix <- model.matrix(f, D$meta_data)

    res <- aldex(t(D$count_data), model_matrix, denom = "all", test="glm")
    
    res <- gather(as.data.frame(res) %>% rownames_to_column('taxa'), "measure", "value", 2:(dim(as.data.frame(res))[2]+1))
    
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
    res$feature <- unlist(lapply(res$measure, clean.feature))
    
    clean.determinant <- function(v){
        gsub("^model.", "", gsub(".Estimate", "", gsub(".Pr...t..", "", gsub(".Pr...t...BH", "", gsub(".t.value", "", gsub(".Std..Error", "", v))))))
    }
    res$determinant <- unlist(lapply(res$measure, clean.determinant))

    
    res.full <- merge(merge(merge(coeff, se, by=c("taxa","determinant")), W, by=c("taxa","determinant")), p.val, by=c("taxa","determinant"))
    res.full$formula <- rep(format(f), dim(res.full)[1])
    res.full$method <- rep("ancombc", dim(res.full)[1])
    res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
    res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

    res <- res.full
    res <- res[startsWith(res$determinant, mainvar),]
    res <- res[res$taxa %in% D$nonrare,]

    res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
    res$determinant <- rep(mainvar, dim(res)[1])
    names(res)[names(res)=="coefficient"] <- "effectsize"

    return(list(res=res, res.full=res.full))
}

R <- lapply(1:length(D$formula$main_var), do)

res <- bind_rows(lapply(R, function(x){x$res}))
res$qvalue <- p.adjust(res$pvalue, "fdr")
res.full <- bind_rows(lapply(R, function(x){x$res.full}))

###############################################################################
# Output

write_tsv(res, D$outfile)
write_tsv(res.full, paste0(c(D$outfile, ".full.tsv"), collapse=""))