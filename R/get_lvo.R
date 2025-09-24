#' Helper function to obtain elements from a fitted lavaan object
#' @keywords internal
#' @param fit see gen_fit2
#' @noRd
get_lvo <- function(fit) {
  pt <- lavaan::parTable(fit)
  if (fit@Options$model.type == "cfa") {
    pt <- pt %>% dplyr::filter(op == "=~")
    n.f <- unique(pt %>% dplyr::select(lhs))$lhs
    mod <- vector("list", length(n.f))
    for (i in 1:length(n.f)) {
      mod[[i]] <-
        paste0(n.f[i], "=~", paste0((pt %>% dplyr::filter(lhs == n.f[i]))$rhs, collapse = "+"))
    }
    mod <- paste0(unlist(mod), collapse = "\n")
  }
  if (fit@Options$model.type != "cfa") {
    n.ov <-
      unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(rhs))$rhs
    n.f <-
      unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(lhs))$lhs
    lam.mod <- vector("list", length(n.f))
    for (i in 1:length(n.f)) {
      lam.mod[[i]] <-
        paste0(n.f[i], "=~", paste0((pt %>% dplyr::filter(op == "=~") %>% dplyr::filter(lhs == n.f[i]))$rhs, collapse = "+"))
    }
    ptt <- pt %>% dplyr::filter(op != "=~", free != 0)
    par.mod <- mapply(function(lhs, op, rhs) {
      paste0(lhs, op, rhs)
    }, ptt$lhs, ptt$op, ptt$rhs)
    mod <- paste0(c(unlist(lam.mod), par.mod), collapse = "\n")
  }
  x <- data.frame(lavaan::lavInspect(fit, "data"))
  n <- lavaan::lavInspect(fit, "ntotal")
  esti <- fit@Options$estimator
  return(res = list(
    x = x,
    n = n,
    mod = mod,
    esti = esti
  ))
}
