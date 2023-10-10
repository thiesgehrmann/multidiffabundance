###############################################################################
# Run CORNCOB

mda.corncob <- function(mda.D, ...){
    D <- mda.D
    
    suppressPackageStartupMessages({
        require(corncob)
        require(phyloseq)
        require(tibble)
        require(dplyr)})

    do <- function(f_idx){
        fdata <- D$formula[[f_idx]]
        f <- fdata$fn
        
        if ( length(fdata$parts.random) > 0 ){
            message(paste0(c("[MDA] mda.corncob: Formula on ", f_idx, " contains random effects. Corncob can not handle random effects")))
            return(mda.common_do(D, f_idx, mda.empty_output(D, f_idx, "Formula incompatible with Corncob analysis (random effect specified)"), "corncob", skip_taxa_sel=TRUE))
        }

        results <- mda.cache_load_or_run_save(D, f_idx, "corncob",{
                       meta_data.nona <- na.omit(fdata$data)
                       count_data.nona <- D$count_data[rownames(meta_data.nona),]
            
                       OTU <- phyloseq::otu_table(t(count_data.nona), taxa_are_rows = T)
                       sampledata <- phyloseq::sample_data(meta_data.nona, errorIfNULL = F)
                       phylo <- phyloseq::merge_phyloseq(OTU, sampledata)
            
                       corncob::differentialTest(formula= f,
                                                 phi.formula = f,
                                                 phi.formula_null = f,
                                                 formula_null = ~ 1,
                                                 test="Wald", data=phylo,
                                                 boot=F,
                                                 fdr_cutoff = 0.05) })

        ram <- results$all_models
        ram <- lapply(1:length(ram), function(t_idx){
            ms <- ram[[t_idx]]
            if (!('coefficients' %in% names(ms))){
                return(NULL)
            }
            ms <- as.data.frame(ms$coefficients)
            ms <- ms[startsWith(rownames(ms), "phi."),]
            ms$taxa <- rep(colnames(D$count_data)[t_idx], dim(ms)[1])
            ms <- ms %>% rownames_to_column('variable.mda')
            ms$variable.mda <- gsub("^phi[.]", "", ms$variable.mda)

            names(ms)[names(ms)=="Estimate"] <- "effectsize"
            names(ms)[names(ms)=="Std. Error"] <- "se"
            names(ms)[names(ms)=="t value"] <- "stat"
            names(ms)[names(ms)=="Pr(>|t|)"] <- "pvalue"
            ms
        })

        res.full <- bind_rows(ram)
        
        mda.common_do(D, f_idx, res.full, "corncob", skip_taxa_sel=TRUE)

    }

    R <- lapply(1:length(D$formula), do)

    mda.common_output(R)
}