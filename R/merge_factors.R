#' Internal function that takes a model and merges factors for discriminant validity testing (merging)
#'
#' @param fit A fitted lavaan model.
#' @param merged.factors Names of the factors to be merged. Must be equal to 2. If missing (the default), the first and second factor of the model are selected.
#' The first factor named will be retained while the second factor will be dropped.
#' @keywords internal
#' @return A lavaan parameter table of the merged factors.
#' @noRd
merge_factors <- function(fit, merged.factors = NULL) {
  if (!inherits(fit, "lavaan"))
    stop("Object fit is not a fitted model from lavaan. Please revise.")
  checkmate::assertVector(merged.factors, len = 2, null.ok = TRUE)
  std.lv <- ifelse(fit@Options$std.lv, TRUE, FALSE)
  pt <- as.data.frame(fit@ParTable)
  # vars <-
  #   unique(pt[which(grepl("^[[:upper:]]", pt$lhs) == TRUE &
  #                     grepl("^[[:upper:]]", pt$rhs) == TRUE), "lhs"])
  vars <- rownames(lavaan::inspect(fit, "veta"))
  if (!is.null(merged.factors))
    mf <- match(merged.factors, vars)
  if (is.null(merged.factors))
    mf <- c(1, 2)
  if (any(is.na(mf)))
    stop("At least one of the factors to be merged is not a factor in the model. Please revise.")
  pt$lhs[which(pt$lhs == vars[mf[2]] &
                 pt$op == "=~")] <- vars[mf[1]]
  pt2 <-
    lapply(pt, function(x)
      x[seq_len(length(pt$lhs))[-which(pt$rhs == vars[mf[2]] |
                                         pt$lhs == vars[mf[2]])]])
  pt2$id <- seq_len(length(pt2$id))
  if (std.lv) {
    fw <- which(pt2$free == 0)
  }
  if (!std.lv) {
    lw <- which(pt2$lhs == vars[mf[1]] & pt2$op == "=~")
    pt2$ustart[lw] <- c(1, rep(NA, length(lw) - 1))
    fw <- which(pt2$free == 0)[-2]
  }
  pt2$free[-fw] <- seq(1:(length(pt2$free) - length(fw)))
  #(seq_len(length(pt2$free)) - length(fw))
  pt2$start <- NULL
  pt2$est <- NULL
  pt2$se <- NULL
  pt2 <- as.data.frame(pt2)
  return(pt2)
}
