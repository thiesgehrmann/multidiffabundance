 ###############################################################################
# Group vs group analysis - CST analysis

mda.group <- function(mda.D, group.col="group", ...){
    D <- mda.D
    
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tibble))
    
    Dcopy <- D
    Dcopy$count_data <- D$count_data[!is.na(D$meta_data[,group.col]),]
    Dcopy$meta_data <- D$meta_data[!is.na(D$meta_data[,group.col]),]
    D <- Dcopy
    
    lmgroup <- function(D, formula){

        f <- update(formula, mda.group.col ~ .)
        
        res <- lapply(unique(D$meta_data[,group.col]), function(g){
            meta_data <- D$meta_data
            meta_data$mda.group.col <- as.numeric(D$meta_data[,group.col] == g)
            fit <- glm(f, data=meta_data, na.action = 'na.exclude', family="binomial")
            s <- as.data.frame(coefficients(summary(fit)))
            s$taxa <- rep(g, dim(s)[1])
            s <- s %>% rownames_to_column("variable")
            s
        })
        res <- bind_rows(res)

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="z value"] <- "stat"
        names(res)[names(res)=="Pr(>|z|)"] <- "pvalue"

        res
    }
    
    lmergroup <- function(D, formula){
        suppressPackageStartupMessages(library(lmerTest))

        f <- update(formula, mda.group.col ~ .)

        res <- lapply(unique(D$meta_data[,group.col]), function(g){
            meta_data <- D$meta_data
            meta_data$mda.group.col <- D$meta_data[,group.col] == g
            fit <- glmer(f, data=meta_data, na.action = 'na.exclude', family="binomial")
            s <- as.data.frame(coefficients(summary(fit)))
            s$taxa <- rep(g, dim(s)[1])
            s <- s %>% rownames_to_column("variable")
            s
        })
        res <- bind_rows(res)

        names(res)[names(res)=="Estimate"] <- "effectsize"
        names(res)[names(res)=="Std. Error"] <- "se"
        names(res)[names(res)=="z value"] <- "stat"
        names(res)[names(res)=="Pr(>|z|)"] <- "pvalue"

        res
    }

    do <- function(f_idx){

        method <- if ( (length(D$formula$rand_intercept[[f_idx]]) + length(D$formula$rand_slope[[f_idx]])) > 0 ){
            lmergroup
        } else { lmgroup }
        f <- D$formula$formula[[f_idx]]
        mainvar <- D$formula$main_var[f_idx]

        res.full <- mda.cache_load_or_run_save(D, "group", f, method(D, f))

        res.full$formula <- rep(mda.deparse(f), dim(res.full)[1])
        res.full$method <- rep("group", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

        # Select only the relevant variable ( taxa are selected already in lmclr, but we repeat it here for safety )
        res <- res.full[startsWith(res.full$variable, mainvar),]

        res$qvalue.withinformula <- p.adjust(res$pvalue, "fdr")
        res$variable <- rep(mainvar, dim(res)[1])

        return(list(res=res, res.full=res.full))
    }

    R <- lapply(1:length(D$formula$main_var), do)


    res <- bind_rows(lapply(R, function(x){x$res}))
    res$qvalue <- p.adjust(res$pvalue, "fdr")
    res.full <- bind_rows(lapply(R, function(x){x$res.full}))

    ###############################################################################
    # Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}