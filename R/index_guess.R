#' Helper function that guesses GoF or BoF from a given index name
#'
#' @param index A fit index or measure provided by function fitmeasures in package lavaan
#' @keywords internal
#' @return Returns GoF (Goodness-of-Fit index) or BoF (Badness of Fit index).
#' @examples
#' index_guess("cfi")
#' index_guess("tli")
#' index_guess("rmsea")
#' index_guess("srmr")
#' @noRd
index_guess <- function(index) {
  #Expand later
  bof <- c("rmsea", "rmr", "srmr", "crmr")
  gof <-
    c("cfi",
      "tli",
      "nnfi",
      "rfi",
      "nfi",
      "pnfi",
      "ifi",
      "rni",
      "gfi",
      "agfi",
      "pgfi",
      "mfi")
  if (grepl(".", index)) {
    index <- strsplit(index, "[.]")[[1]][1]
  }
  if (grepl("_", index)) {
    index <- strsplit(index, "[_]")[[1]][1]
  }
  idx <- tolower(index)
  r <-
    ifelse(idx %in% gof, "GoF", ifelse(idx %in% bof, "BoF", "not a fit index"))
  return(r)
}
