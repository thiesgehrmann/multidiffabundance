###############################################################################
# Maaslin2

mda.maaslin2 <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages(require(Maaslin2))
    suppressPackageStartupMessages(require(dplyr))
    suppressPackageStartupMessages(require(tibble))

    do <- function(f_idx){
        f <- D$formula$norand[[f_idx]]
        f.cache <- f
        mainvar <- D$formula$main_var[f_idx]
        
        if ( length(D$formula$rand_slope[[f_idx]]) > 0 ){
            message(paste0(c("[MDA] mda.maaslin2: Formula on ", f_idx, " contains random slope effects. Maaslin2 can not handle random slopes. Run will continue without random slopes")))
        }
        
        random_effects <- NULL
        if ( length(D$formula$rand_intercept[[f_idx]]) > 0 ){
            randvars <- D$formula$rand_intercept[[f_idx]]
            random_effects <- unlist(lapply(randvars, function(x){ trimws(gsub(")", "", unlist(strsplit(x, split="\\|"))[[2]]))}))
            f.cache <- update.formula(f.cache, paste0(c("~.+", paste0(D$formula$rand_intercept[[f_idx]], collapse="+")), collapse=""))
        }
        
        fix_factors <- D$meta_data
        possible_factors <- colnames(fix_factors)[unlist(lapply(colnames(fix_factors), function(x) { typeof(fix_factors[,x])})) == 'character']
        dumb_masslin_factor_crap <- function(d){
            fd <- factor(fix_factors[,d])
            lv <- sort(levels(fd))
            if (length(lv) > 2){
                paste0(c(d, lv[1]), collapse=",")
            } else {
                NULL
            }
        }
        maaslin2.reference <- paste0(unlist(lapply(possible_factors, dumb_masslin_factor_crap)), collapse=";")

        mas <- mda.cache_load_or_run_save(D, "maaslin2", f.cache,
                    Maaslin2(input_data = D$count_data,
                             input_metadata = D$meta_data,
                             output = paste0(c(D$outprefix, "/maaslin2.output.folder"), collapse=""),
                             min_abundance = 0.0,
                             min_prevalence = 0.0,
                             normalization = "TSS",
                             transform = "LOG",
                             analysis_method = "LM",
                             max_significance = 0.05,
                             fixed_effects = labels(terms(f)),
                             random_effects = random_effects,
                             correction = "BH",
                             standardize = FALSE,
                             plot_heatmap = FALSE,
                             plot_scatter = FALSE,
                             cores = 1,
                             reference = maaslin2.reference))

        res.full <- as_tibble(as.data.frame(mas$results))
        res.full <- res.full[,c("feature","metadata","coef","stderr","pval")]
        res.full$formula <- rep(mda.deparse(f.cache), dim(res.full)[1])
        res.full$method <- rep("maaslin2", dim(res.full)[1])
        res.full$n <- rep(mda.meta.n(D, mainvar), dim(res.full)[1])
        res.full$freq <- rep(mda.meta.freq(D, mainvar), dim(res.full)[1])
        res.full$stat <- rep(NA, dim(res.full)[1])

        names(res.full)[names(res.full)=="feature"] <- "taxa"
        names(res.full)[names(res.full)=="metadata"] <- "variable"
        names(res.full)[names(res.full)=="coef"] <- "effectsize"
        names(res.full)[names(res.full)=="stderr"] <- "se"
        names(res.full)[names(res.full)=="pval"] <- "pvalue"

        res <- res.full[(res.full$variable == mainvar),]
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
    # Output

    column.order <- c("taxa","variable","effectsize","se","stat","pvalue","qvalue.withinformula","qvalue","formula","method","n","freq")
    res <- res[,column.order]
    return(list(res=res, res.full=res.full, summary=mda.summary(res)))
}