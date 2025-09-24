#' Internal function that checks standardization and changes cutoff if needed
#'
#' @param fit A fitted lavaan model.
#' @param dv.factors ...
#' @param dv.cutoff ...
#' @return A cutoff for discriminant validity testing if model is not standardized
#' @keywords internal
#' @noRd
check_std <- function(fit, dv.factors, dv.cutoff) {
  psi <- lavaan::inspect(fit, "free")$psi
  std <- any(diag(psi) == 0)
  if (!std) {
    sc <-
      (lavaan::inspect(fit, "cor.lv") / lavaan::inspect(fit, "cov.lv"))[dv.factors[1], dv.factors[2]]
    dv.cutoff <- dv.cutoff / sc
  }
  return(round(dv.cutoff, 3))
}
