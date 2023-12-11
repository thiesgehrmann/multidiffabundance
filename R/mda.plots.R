#' Return an interactive differential abundance heatplot of the effect sizes of the covariates on the taxa using different methods
#'
#' @param res The result of an mda.all(D, alpha=TRUE) or mda.alpha run.
#' @param low The color to use for low effect sizes.
#' @param high The color to use for high effect sizes.
#' @param mid The color to use for effect sizes around 0.
#' @export
mda.plot.heatplot <- function(res, ta=NULL, low="magenta", high="limegreen", mid="white", min_sign=1, min_eff, q_lim=0.05) {

  if (!is.null(ta)) {
    require("tidytacos")
    # Rename taxa using tidytacos obj
    res$res.full <- res$res.full %>% 
      dplyr::left_join(
        ta %>% 
        tidytacos::add_taxon_name() %>% 
        tidytacos::taxa() %>% 
        dplyr::select(c(taxon_name, taxon_id)), 
        by=c("taxa"="taxon_id")) %>% 
      dplyr::select(-taxa) %>%
      dplyr::rename(taxa=taxon_name)
  }

  qhits <- res$res.full %>% 
    dplyr::filter(qvalue<q_lim, !method %in% c("alpha", "beta")) %>% 
    dplyr::filter(!is.na(variable), variable != "(Intercept)") %>% 
    dplyr::group_by(method, variable, taxa) %>% 
    dplyr::summarize(max(qvalue)) %>% 
    dplyr::group_by(taxa, variable) %>% 
    dplyr::summarize(nq=dplyr::n())

  res$res.full <- res$res.full %>% 
    dplyr::left_join(qhits)
  
  res$res.full$nq[is.na(res$res.full$nq)] <- 0

  sign_tax <- unique(res$res.full %>% 
  dplyr::filter(
    !is.na(variable), 
    variable != "(Intercept)") %>%
  dplyr::filter(nq >= min_sign) %>% 
  dplyr::filter(!is.na(qvalue)) %>%
  dplyr::filter(abs(effectsize) >= min_eff) %>%
  dplyr::pull(taxa))
  
  if (length(sign_tax) == 0) {
    stop("No taxa meet the criteria for significance.")
  }

  sign_res <- res$res.full %>% 
  dplyr::filter(taxa %in% sign_tax) %>%
  dplyr::filter(
    !is.na(variable), 
    variable != "(Intercept)",
    !method %in% c("alpha", "beta"))

  methods <- unique(sign_res$method)
  make_button <- function(method) {
    list(
      method = "restyle",
      args = list("transforms[0].value", method),
      label = method
    )
  }
  buttons <- lapply(methods, make_button)

  plotly::plot_ly(
    sign_res,
    x=~taxa,
    y=~variable,
    z=~effectsize,
    zmid=0,
    type="heatmap",
    colors = colorRamp(c(low,mid,high)),
    text=~qvalue,
    transforms = list(
      list(
        type = 'filter',
        target = ~method,
        operation = '==',
        value = "maaslin2"
      )
    )
  ) %>% plotly::add_annotations(
    x=sign_res$taxon_name,
    y=sign_res$variable,
    text=ifelse(
      is.na(sign_res$nq), "", sign_res$nq), 
      showarrow=FALSE) %>% plotly::layout(
    updatemenus = list(
      list(
        yanchor= "bottom",
        direction= "right",
        type = 'dropdown',
        active = 1,
        buttons = buttons
      )
    )
  )

}

#' Return an interactive plot of the effect sizes of the covariates on the beta diversity
#'
#' @param res The result of an mda.all(D, alpha=TRUE) or mda.alpha run.
#' @export
mda.plot.beta <-  function(res) {
  res.beta <- res$res.full %>%
    dplyr::filter(
      method == "beta",
      !is.na(variable))

  plotly::plot_ly(
    res.beta,
    x = ~effectsize,
    y = ~variable,
    type = "bar",
    text = ~paste("P-value:", round(pvalue,4))
  )
}

#' Return an interactive plot of the effect sizes of the covariates on the alpha diversity
#'
#' @param res The result of an mda.all(D, alpha=TRUE) or mda.alpha run.
#' @export
mda.plot.alpha <-  function(res) {
  res.alpha <- res$res.full %>% 
  dplyr::filter(
    method == "alpha", 
    variable != "(Intercept)",
    !is.na(variable)
  )
  
  plotly::plot_ly(
  res.alpha, 
  x=~effectsize, 
  y=~variable, 
  type="scatter", mode="markers",
  error_x=~list(array=se, color="#000000"),
  text=~paste("P-value:", round(pvalue,4))
)
}

#' Return a heatplot, beta plot, and alpha plot in a single plotly object.
#'
#' @param res The result of an mda.all(D, alpha=TRUE) or mda.alpha run.
#' @param low The color to use for low effect sizes.
#' @param high The color to use for high effect sizes.
#' @param mid The color to use for effect sizes around 0.
#' @export
mda.plot.all <- function(res, ...) {

  hp <- mda.plot.heatplot(res, ...)
  bp <- mda.plot.beta(res)
  ap <- mda.plot.alpha(res)

  plotly::subplot(
    hp, ap, bp, 
    nrows=1, shareY = TRUE
  ) %>% 
  plotly::layout(
    annotations = list(
       list(x=0.2, y=1, text="Alpha diversity", showarrow=F, xref="x2", yref="paper", font= list(size=18)),
       list(x=0.2, y=1, text="Beta diversity", showarrow=F, xref="x3", yref="paper", font = list(size=18))
    ), width = 1000, height = 600
  ) %>%
  plotly::config(
  modeBarButtonsToRemove = c("zoomIn2d", "zoomOut2d", "zoom2d", "pan2d", "lasso2d", "select2d")
  )
}
