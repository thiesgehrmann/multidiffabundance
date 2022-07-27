###############################################################################
# Run BETA

mda.beta <- function(mda.D, permutations=999){
    require(vegan)
    require(tidyr)
    require(tibble)
    
    D <- mda.D

    do <- function(f_idx){
        f <- D$formula$norand[[f_idx]]
        mainvar <- D$formula$main_var[f_idx]
        
        if ( length(D$formula$rand_slope[[f_idx]]) > 0 ){
            message(paste0(c("[MDA] mda.beta. Formula on ", f_idx, " contains random slope effects. Adonis2 can not handle random slopes. Run will continue without random slopes")))
        }
        
        block <- NULL
        if ( length(D$formula$rand_intercept[[f_idx]]) > 0 ){
            rblock <- D$formula$rand_intercept[[f_idx]][1]
            block <- trimws(gsub(")", "", unlist(strsplit(rblock, split="\\|"))[[2]]))
            message(paste0(c("[MDA] mda.beta: Formula on ", f_idx, " contains a random intercept effect. Adonis2 cannot handle random effects. However, we will translate this into a block design. This means that permutations will be performed within these blocks. If you do not desire this, replace the random intercept with a fixed effect. Will block with (1|", block, ").")))
            if (length(D$formula$rand_intercept[[f_idx]]) > 1){
                message(paste0(c("[MDA] mda.beta: Formula on ", f_idx, " contains more than one random intercept effect. Adonis2 can only handle one random effect. Will continue with (1|", block, ").")))
            }
        }
        f.cache <- if (is.null(block)){f}else{paste0(c(mda.deparse(f), " + (1|", block, ")"), collapse="")}
        f.cache <- as.formula(f.cache)
        variables <- c(mainvar, D$formula$adj_var[[f_idx]], block)

        meta_data.nona <- na.omit(D$meta_data[,variables,drop=FALSE])
        count_data.nona <- D$count_data[rownames(meta_data.nona),]

        Dcopy <- D
        Dcopy$count_data <- count_data.nona
        
        out <- mda.cache_load_or_run_save(Dcopy, "beta", f.cache, {
            meta_data.nona <- na.omit(D$meta_data[,variables,drop=FALSE])
            count_data.nona <- D$count_data[rownames(meta_data.nona),]
            
            dist <- eval(mda.cache_load_or_run_save(Dcopy, "beta_dist", f.cache, {vegdist(count_data.nona)}))
            strata <- if (is.null(block)){NULL}else{meta_data.nona[,block]}
            f <- update(f, dist ~ .)
            assign("dist",dist,envir=environment(f)) # Dumbest shit I ever saw
            
            adonis2(f, meta_data.nona, na.action='na.exclude', strata=strata, by="margin", permutations=permutations)
        })

        res.full <- out
        res.full <- res.full[,c("R2","F", "Pr(>F)")]


        names(res.full)[names(res.full)=="R2"] <- "effectsize"
        names(res.full)[names(res.full)=="Std. Error"] <- "se"
        names(res.full)[names(res.full)=="F"] <- "stat"
        names(res.full)[names(res.full)=="Pr(>F)"] <- "pvalue"
        
        print(dim(res.full)[1])
        res.full$variable <- rownames(res.full)
        res.full$se <- rep(NA, dim(res.full)[1])
        res.full$taxa <- rep("mda.beta", dim(res.full)[1])
        res.full$qvalue.withinformula <- res.full$pvalue
        
        res.full$formula <- rep(mda.deparse(f), dim(res.full)[1])
        res.full$method <- rep("beta", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])

        res <- res.full[startsWith(res.full$variable, mainvar),]
        res$variable <- rep(mainvar, dim(res)[1])

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