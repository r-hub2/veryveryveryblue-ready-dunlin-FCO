#' Helper function that gets free model from a population model.
#' @param pop.mod A population model.
#' @keywords internal
#' @return The syntax of a free lavaan model
#' @noRd
get_free_mod <- function(pop.mod) {
  pt <- lavaan::lavaanify(pop.mod)
  pt <- pt %>% dplyr::filter(op == "=~")
  n.f <- unique(pt %>% dplyr::select(lhs))$lhs
  free <- vector("list", length(n.f))
  for (i in 1:length(n.f)) {
    free[[i]] <-
      paste0(n.f[i], "=~", paste0((pt %>% dplyr::filter(lhs == n.f[i]))$rhs, collapse = "+"))
  }
  free <- paste0(unlist(free), collapse = "\n")
  return(free)
}
