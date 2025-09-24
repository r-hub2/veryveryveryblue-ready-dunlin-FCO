#' Helper function to get power and related stats
#'
#' @keywords internal
#' @param cutoff the cutoff to be used
#' @param index fit index
#' @param correct.fits see flex_co2
#' @param miss.fits ...
#' @param null.is.miss assuming none or missing values to be misfit?
#' @noRd
power.finder <-
  function(cutoff,
           index,
           correct.fits,
           miss.fits,
           NULL.is.miss = T) {
    fi <- as.character(index)
    cutoff <- as.numeric(cutoff)
    ix <- index_guess(toupper(fi))
    cf <- correct.fits
    mf <- miss.fits
    #True negative: Correct models correctly confirmed by cutoff (healthy)
    #False positive: Correct models wrongly rejected by cutoff (ill)
    #True positive: Misspec. models correctly rejected by cutoff (ill)
    #False negative: Misspec. models wrongly confirmed by cutoff (healthy)
    if (ix == "GoF") {
      if (NULL.is.miss) {
        cx <- dplyr::if_else(is.na(cutoff), 1, cutoff)
      }
      if (!NULL.is.miss) {
        cx <- dplyr::if_else(is.na(cutoff), 0, cutoff)
      }
      neg <-
        dplyr::if_else(dplyr::pull(cf[, tolower(fi)]) >= cx, "TN", "FP")
      pos <-
        dplyr::if_else(dplyr::pull(mf[, tolower(fi)]) < cx, "TP", "FN")
    }
    if (ix == "BoF")
    {
      if (NULL.is.miss) {
        cx <- dplyr::if_else(is.na(cutoff), 0, cutoff)
      }
      if (!NULL.is.miss) {
        cx <- dplyr::if_else(is.na(cutoff), 1, cutoff)
      }
      neg <-
        dplyr::if_else(dplyr::pull(cf[, tolower(fi)]) <= cx, "TN", "FP")
      pos <-
        dplyr::if_else(dplyr::pull(mf[, tolower(fi)]) > cx, "TP", "FN")
    }
    tn <- length(neg[neg == "TN"])
    fp <- length(neg[neg == "FP"])
    tp <- length(pos[pos == "TP"])
    fn <- length(pos[pos == "FN"])
    type1 <- fp / (fp + tn)
    type2 <- fn / (tp + fn)
    type_sum <- type1 + type2
    power <- tp / (tp + fn)
    spec <- tn / (fp + tn)
    res <- c(tn, fp, tp, fn, type1, type2, type_sum, power, spec)
    names(res) <-
      c("TN",
        "FP",
        "TP",
        "FN",
        "TypeI",
        "TypeII",
        "SumTypes",
        "Power",
        "Specificity")
    return(res)
  }
