###############################################################################
#' Maaslin2
#' @export
mda.maaslin2 <- function(mda.D, ...){
    D <- mda.D
    suppressPackageStartupMessages({
        require(Maaslin2)
        require(dplyr)
        require(tibble)})

    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn
        
        if ( length(fdata$parts.random.slope) > 0 ){
            message(paste0(c("[MDA] mda.maaslin2: Formula ", f_idx, " contains random slope effects. Maaslin2 can not handle random slopes.")))
            return(mda.common_do(D, f_idx, mda.empty_output(D, f_idx, "Formula incompatible with maaslin2 analysis (random slope specified)"), "maaslin2", skip_taxa_sel=TRUE))
        }
        
        random_effects <- NULL
        if ( length(fdata$parts.random.intercept) > 0 ){
            random_effects <- unlist(lapply(strsplit(fdata$parts.random.intercept, split="\\|"), function(v){trimws(v[2])}))
        }
        
        fix_factors <- fdata$data
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
        
        nn.metadata <- tidyr::drop_na(fdata$data)
        nn.counts   <- D$count_data[rownames(nn.metadata),]
        
        maaslin2.reference <- paste0(unlist(lapply(possible_factors, dumb_masslin_factor_crap)), collapse=";")
        
        mas <- mda.cache_load_or_run_save(D, f_idx, "maaslin2",
            rlang::exec(
                Maaslin2,
                input_data = nn.counts,
                input_metadata = nn.metadata,
                output = paste0(c(D$outprefix, "/maaslin2.output.folder"), collapse=""),
                fixed_effects = fdata$parts.fixed,
                random_effects = random_effects,
                reference = maaslin2.reference,
                !!!D$args$maaslin2
            ))
        


        res.full <- as_tibble(as.data.frame(mas$results))
        res.full <- res.full[,c("feature","metadata","coef","stderr","pval")]
        names(res.full)[names(res.full)=="feature"] <- "taxa"
        names(res.full)[names(res.full)=="metadata"] <- "variable.mda"
        names(res.full)[names(res.full)=="coef"] <- "effectsize"
        names(res.full)[names(res.full)=="stderr"] <- "se"
        names(res.full)[names(res.full)=="pval"] <- "pvalue"
        
        res.full <- as.data.frame(res.full)
        
        res <- mda.common_do(D, f_idx, res.full, "maaslin2", skip_taxa_sel=FALSE)
        
        res$res.full$se <- res$res.full$se / sqrt(as.numeric(res$res.full$n))
        res$res$se <- res$res$se / sqrt(as.numeric(res$res$n))
        res
    }

    R <- lapply(1:length(D$formula), do)

    mda.common_output(R)
}
