###############################################################################
# ALDEx2

mda.aldex2 <- function(mda.D){
    require("ALDEx2")
    require(tidyr)
    require(tibble)
    D <- mda.D

    do <- function(f_idx){
        f <- D$formula$norand[[f_idx]]
        mainvar <- D$formula$main_var[f_idx]

        if ( (length(D$formula$rand_intercept[[f_idx]]) + length(D$formula$rand_slope[[f_idx]])) > 0 ){
            message(paste0(c("[MDA] mda.aldex2: Formula ", f_idx, " contains random effects. ALDEx2 can not handle random effects. Run will continue without random effects.")))
        }

        out <- mda.cache_load_or_run_save(D, "aldex2", f, {
            mm <- model.matrix(f, D$meta_data)
            counts <- t(D$count_data[rownames(mm),]) # Remove the NAs that have been removed by model.matrix
            aldex(counts, mm, denom = "all", test="glm")
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

        clean.variable <- function(v){
            gsub("^model.", "", gsub(".Estimate", "", gsub(".Pr...t..", "", gsub(".Pr...t...BH", "", gsub(".t.value", "", gsub(".Std..Error", "", v))))))
        }
        res.full$variable <- unlist(lapply(res.full$measure, clean.variable))

        res.full <- pivot_wider(res.full, id_cols=c("taxa", "variable"), names_from=feature, values_from=value)
        res.full$formula <- rep(mda.deparse(f), dim(res.full)[1])
        res.full$method <- rep("aldex2", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])


        res <- res.full[startsWith(res.full$variable, mainvar),]
        res$variable <- rep(mainvar, dim(res)[1])
        res <- res[res$taxa %in% D$nonrare,]

        res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")

        return(list(res=res, res.full=res.full))
    }

    R <- lapply(1:length(D$formula$main_var), do)

    res <- bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    res.full <- bind_rows(lapply(R, function(x){x$res.full}))

    ###############################################################################
    # Format Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}