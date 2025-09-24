#' Internal new, faster generator for gen_fit2
#' @keywords internal
#' @param pop.mod see gen_fit2
#' @param seed ...
#' @param mod ...
#' @param n ...
#' @param sk ...
#' @param ku ...
#' @param esti ...
#' @noRd
# mini.generator <- function(pop.mod, seed, mod, n, sk, ku, esti = "ML") {
#   y <- try(lavaan::simulateData(
#     model = pop.mod,
#     model.type = "cfa",
#     sample.nobs = n,
#     skewness = sk,
#     kurtosis = ku,
#     seed = seed,
#     auto.fix.first = FALSE,
#     std.lv = TRUE
#   ),
#   silent = TRUE)
#   sf <- NA
#   rf <- rep(NA, 46)
#   if (is.data.frame(y)) {
#     sf <- try(lavaan::cfa(
#       mod,
#       y,
#       estimator = esti,
#       auto.fix.first = FALSE,
#       std.lv = TRUE,
#       warn = FALSE
#     ),
#     silent = TRUE)
#   }
#   if (typeof(sf) == "S4") {
#     if (sf@Fit@converged) {
#       rf <- lavaan::fitmeasures(sf)
#     }
#   }
#   return(rf)
# }
# mini.generator <- compiler::cmpfun(mini.generator)
mini.generator <- function(pop.mod, seed, mod, n, sk, ku, esti = "ML") {
  rf <- rep(NA, 46)  # leeres Ergebnis
  y <- try(lavaan::simulateData(
    model = pop.mod,
    model.type = "cfa",
    sample.nobs = n,
    skewness = sk,
    kurtosis = ku,
    seed = seed,
    auto.fix.first = FALSE,
    std.lv = TRUE
  ), silent = TRUE)
  if (inherits(y, "try-error")) return(rf)
  sf <- try(lavaan::cfa(
    mod,
    y,
    estimator = esti,
    auto.fix.first = FALSE,
    std.lv = TRUE,
    warn = FALSE
  ), silent = TRUE)
  if (inherits(sf, "try-error")) return(rf)
  if (sf@Fit@converged) {
    rf <- lavaan::fitmeasures(sf)
  }
  return(rf)
}
