#' Write an mda object to two tsv tables and a file containing the formulas in a folder
#'
#' @param mda the mda object
#' @param folder the folder to save the tables and formulas to
#' @export
mda.write <- function(mda, folder='mda') {
  
  count_dataf <- paste0(folder, "/countdata.tsv")
  meta_dataf <- paste0(folder, "/metadata.tsv")
  formulaf <- paste0(folder, "/formulas.txt")

  fmls <- trimws(paste0(purrr::map(mda$formula, "fn.raw"), collapse="\n"))
  counts <- data.frame(mda$count_data) %>% tibble::rownames_to_column("sample")

  dir.create(folder, showWarnings=FALSE)
  write.table(counts, file=count_dataf, quote=FALSE, sep="\t", row.names=FALSE)
  write.table(mda$meta_data, file=meta_dataf, quote=FALSE, sep="\t", row.names=FALSE)
  cat(fmls, file=formulaf)

}
