 ###############################################################################
# Wilcoxon
mda.wilcoxon <- function(mda.D, wilcoxon.norm="clr", wilcoxon.resid=TRUE, ...){
    D <- mda.D
    
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tibble))
    library(lme4)

    norm_data <- if (wilcoxon.norm == 'clr'){
        as.data.frame(scale(mda.clr(mda.relative_abundance(mda.pseudocount(D$count_data)))))
    } else if (wilcoxon.norm == 'log'){
        as.data.frame(scale(log(mda.relative_abundance(mda.pseudocount(D$count_data)))))
    } else {
        mda.message("The parameter `wilxocon.norm` should be either 'clr' or 'log'.", type="error")
    }

    resids <- function(norm_data, meta_data, formula){
        f <- update(formula, normtaxa ~ .)

        norm_data.resids <- as.data.frame(norm_data)
        meta_data <- as.data.frame(meta_data)

        method <- if ( formula.ismixed(formula) ){ lmer } else { lm }

        norm_data.resids <- apply(norm_data, 2, function(normtaxa){
            meta_data$normtaxa <- normtaxa
            fit <- method(f, data=meta_data, na.action = 'na.exclude')
            if ( formula.ismixed(formula) ){ resid(fit) } else { fit$residuals }
        })
        return(norm_data.resids)
    }
    
    wilcoxon <- function(norm_data, meta_data, formula, taxa=NULL){

        taxa <- if (is.null(taxa)) colnames(count_data) else taxa

        
        var1 <- formula.parts(formula)[1]

        input_data <- if ((length(formula.parts(formula)) > 1) & wilcoxon.resid){
            fn.regress_out <- update(formula, as.formula(paste0(c(" ~ . - ", formula.parts(formula)[1]), collapse="")))

            resids(norm_data[,taxa], meta_data, fn.regress_out)
            
        } else {
                norm_data
        }


        res <- lapply(taxa, function(t){
            g <- meta_data[,var1]
            g1.val <- min(unique(g))
            g2.val <- max(unique(g))

            wt <- wilcox.test(input_data[g == g1.val,t], input_data[g == g2.val,t])

            nvars <- length(formula.parts.fixed(formula))

            data.frame(variable.mda=formula.parts.fixed(formula),
                       effectsize=NA,
                       se=NA,
                       stat=c(wt$statistic, rep(NA, nvars-1)),
                       pvalue=c(wt$p.value, rep(NA, nvars-1)),
                       comment=c(NA, rep("Not tested by mda.wilcoxon", nvars-1)),
                       taxa=t)
        })
        res <- dplyr::bind_rows(res)

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="t value"] <- "stat"
        names(res)[names(res)=="Pr(>|t|)"] <- "pvalue"

        res
    }

    #return(norm_data)
    
    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]

        var1 <- formula.parts(fdata$fn)[1]

        res.full <- if (length(unique(fdata$data[,var1])) != 2){
            mda.message("mda.wilcoxon: This variable is not binary.", type="error")
            mda.empty_output(D, f_idx, comment="This variable is not binary", taxa=D$nonrare)
        } else {
            mda.cache_load_or_run_save(D, f_idx, "wilcoxon", wilcoxon(norm_data, fdata$data, fdata$fn, D$nonrare), order_invariant=FALSE)
        }
        
        res.full$formula <- rep(mda.deparse(fdata$fn.orig), dim(res.full)[1])
        res.full$method <- rep("wilcoxon", dim(res.full)[1])
        res.full <- left_join(res.full, fdata$nfreq, by="variable.mda")

        # taxa are selected already in wilcoxon, but we repeat it here for safety
        res.full <- res.full[res.full$taxa %in% D$nonrare,]

        # Select only the relevant variable 
        first_var <- formula.parts(fdata$fn.orig)[1]
        res <- res.full[res.full$variable == first_var,]
        res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
        
        res.full$qvalue.withinformula <- p.adjust(res.full$pvalue, "fdr")

        return(list(res=res, res.full=res.full))
    }

    R <- lapply(1:length(D$formula), do)


    res <- dplyr::bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    
    res.full <- dplyr::bind_rows(lapply(R, function(x){x$res.full}))
    res.full$qvalue <- p.adjust(res.full$pvalue, "fdr")

    ###############################################################################
    # Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq","comment")
    res <- res[,column.order]
    res.full <- res.full[, column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}