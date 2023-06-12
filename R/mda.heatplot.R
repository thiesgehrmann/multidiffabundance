mda.heatmap <- function(res, pipeline="maaslin2", ta=NULL, minEff=1, skip_sign=F, low="blue", high="red") {
 
  res_matrix <- function(res, cval="effectsize") {
    #cval <- enquo(col_value)
    res %>% 
    pivot_wider(id_cols=taxa, names_from = variable, values_from = cval) %>%
    column_to_rownames(var="taxa")
  }
  
  # Visualize effectsize/padj for selected method
  if (all(!pipeline %in% res$summary$method_order)) {
    stop(paste("Method", method, "not found in results object."))
  }
  res_df <- res$res %>% dplyr::filter(method == pipeline)
  
  # Filter on effsize
  effsizes_taxa <- res_matrix(res_df) %>% 
    mutate_all(abs) %>%
    mutate(maxEff = do.call(pmax, .) ) %>%
    filter(maxEff >= minEff) %>%
    row.names()
  res_df_min_eff <- res_df %>% filter(taxa %in% effsizes_taxa)
  
  # Filter on significance
  res_signif <- res$summary$nsig %>% column_to_rownames(var="taxa")
  signif_taxa <- res_signif %>% mutate(n_sig = rowSums(.)) %>% 
    filter(n_sig > 0) %>% 
    row.names()
  if (length(signif_taxa) > 0 && !skip_sign) {
  res_df_sign <- res_df_min_eff %>% filter(taxa %in% signif_taxa)
  } 
  else {
    res_df_sign <- res_df_min_eff
    warning("No signifficant results after multiple testing")
  }
  # Optional translate taxon_id to taxnames
  taxids <- res_df_sign$taxa
  if (!is.null(ta)) {
    if (!"taxon_name" %in% colnames(ta$taxa)){
      stop("Tidyamplicons object does not have a taxon_name column.")
    }
    res_df_sign <- res_df_sign %>% 
      left_join(ta$taxa, by=c("taxa" = "taxon_id")) %>% 
      mutate(taxa = taxon_name)
  }
  
  effM <- res_matrix(res_df_sign)
  signM <- res_signif[taxids,]
  signM <- apply(signM, MARGIN = c(1,2), FUN = function(x) {
    ifelse(x == 0, "", x)
  })
  padjM <- res_matrix(res_df_sign, "qvalue")
  padjM <- apply(padjM, MARGIN = c(1,2), FUN = function(x) {
    paste0("Adj P: ", format(signif(x, 3), scientific=T))
  })
  
  heatmaply::heatmaply(
    effM,
    scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(
        low = low, 
        high = high, 
        midpoint = 0),
        show_dendrogram = c(F,F),
        cellnote = signM,
        cellnote_color = "black",
        custom_hovertext = padjM
  )
}
