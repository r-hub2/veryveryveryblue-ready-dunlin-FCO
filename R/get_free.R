#' Helper function that gets free model from a lavaan parameter table.
#'
#' @param pt A lavaan parameter table. Please use factors / latent variables with uppercase names (e.g., F1). The function is an addition to simstandard:fixed2free which only works with models.
#' @param dv.factors The selected factors relevant in case of retaining constraints
#' @param mode Mode to be a relevant for constraints. Only if mode is constraining, constraints are kept.
#' @keywords internal
#' @return The syntax of a free lavaan model
#' @noRd
get_free <- function(pt, dv.factors, mode) {
  # vars <-
  #   unique(pt[which(grepl("^[[:upper:]]", pt$lhs) == TRUE &
  #                     grepl("^[[:upper:]]", pt$rhs) == TRUE), "lhs"])
  vars <-
    unique(pt[which(grepl("^[[:upper:]]", pt$lhs) == TRUE &
                      pt$op == "=~"), "lhs"])
  if (mode == "constraining")
    mod <- rep(NA, length = length(vars) + 1)
  if (mode != "constraining")
    mod <- rep(NA, length = length(vars))
  for (i in seq_len(length(vars))) {
    wh <- which(grepl(vars[i], pt$lhs) == TRUE & pt$op == "=~")
    mod[i] <-
      paste0(vars[i], " =~ ", paste(pt[wh, "rhs"], collapse = " + "))
  }
  if (mode == "constraining") {
    cv <-
      which(
        grepl(dv.factors[1], pt$lhs) == TRUE &
          pt$op == "~~" & grepl(dv.factors[2], pt$rhs) == TRUE
      )
    #Maybe, the dv.factors are reversed
    if (length(cv) == 0) {
      dv.factors <- dv.factors[c(2, 1)]
      cv <-
        which(
          grepl(dv.factors[1], pt$lhs) == TRUE &
            pt$op == "~~" & grepl(dv.factors[2], pt$rhs) == TRUE
        )
    }
    if (length(cv) == 0)
      stop("Constrained correlation for selected factors not found. Perhaps misspelled?")
    mod[length(mod)] <-
      paste0(dv.factors[1], " ~~ ", paste(pt[cv, "ustart"]), "*", dv.factors[2])
  }
  mod <- paste(mod, collapse = "\n")
  return(mod)
}
