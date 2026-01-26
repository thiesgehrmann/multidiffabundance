#' @export
mda.aldex2 <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages({
        require("ALDEx2")
        require(tidyr)
        require(dplyr)
        require(tibble)
        })

    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn

        if ( length(fdata$parts.random) > 0 ){
            message(paste0(c("[MDA] mda.aldex2: Formula on ", f_idx, " contains random effects. Aldex2 can not handle random effects")))
            return(mda.common_do(D, f_idx, mda.empty_output(D, f_idx, "Formula incompatible with Aldex2 analysis (random effect specified)"), "aldex2", skip_taxa_sel=TRUE))
        }

        out <- mda.cache_load_or_run_save(D, f_idx, "aldex2", {
            mm <- model.matrix(f, fdata$data)
            counts <- t(D$count_data[rownames(mm),]) # Remove the NAs that have been removed by model.matrix
            ALDEx2::aldex(counts, mm, denom = "all", test="glm")
        })

        res.full <- tidyr::gather(
            as.data.frame(out) %>% 
            rownames_to_column('taxa'), 
            "measure", "value", 2:(dim(as.data.frame(out))[2]+1))
        r <<- res.full
        clean.feature.1 <- function(v){
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
        
        clean.feature.2 <- function(v){
            ret <- if (endsWith(v, ".Est")) {
                        "effectsize"
                    } else if (endsWith(v, ".SE")) {
                        "se"
                    } else if (endsWith(v, ".t.val")) {
                        "stat"
                    } else if (endsWith(v, ".pval")) {
                        "pvalue"
                    } else if (endsWith(v, ".pval.holm")) {
                        "qvalue"
                    } else if (endsWith(v, ".pval.padj")) {
                        "qvalue"
                    }
            ret
        }
        clean.feature <- if (
            package_version(installed.packages()["ALDEx2", "Version"]) >= package_version("1.3")) {
                clean.feature.2
        } else {clean.feature.1}
        res.full$feature <- unlist(lapply(res.full$measure, clean.feature))
        

        clean.variable.1 <- function(v){
            gsub("^model.", "", 
            gsub(".Estimate", "", 
            gsub(".Pr...t..", "", 
            gsub(".Pr...t...BH", "", 
            gsub(".t.value", "", 
            gsub(".Std..Error", "", v))))))
        }
        
        clean.variable.2 <- function(v){
            gsub("^model.", "", 
            gsub(".Est", "", 
            gsub(".pval", "", 
            gsub(".pval.holm", "", 
            gsub(".t.val", "", 
            gsub(".SE", "", v))))))
        }
        clean.variable <- if (
            package_version(installed.packages()["ALDEx2", "Version"]) >= package_version("1.3")) {
                clean.variable.2
        } else {clean.variable.1}
        res.full$variable.mda <- unlist(lapply(res.full$measure, clean.variable))
        
        res.full <- tidyr::pivot_wider(res.full, id_cols=c("taxa", "variable.mda"), names_from=feature, values_from=value)
        
        mda.common_do(D, f_idx, res.full, "aldex2", skip_taxa_sel=FALSE)

    }

    R <- lapply(1:length(D$formula), do)

    mda.common_output(R)
}
