#' Internal function to generate fits
#' @param x see gen_fit
#' @param n ...
#' @param seed ...
#' @param mode ...
#' @param pop.mod1 ...
#' @param free1 ...
#' @param free2 ...
#' @param s ...
#' @param k ...
#' @param nf ...
#' @keywords internal
#'
#' Allows for one or two models. In case of the latter, it only generates data from the population model of model one.
#' @noRd
generator <- function(x = NULL,
                      n = NULL,
                      seed,
                      mode,
                      pop.mod1,
                      free1 = NULL,
                      free2 = NULL,
                      s = 1,
                      k = 1,
                      nf) {
  y <- NA
  if (is.null(x) & is.null(n)) stop("Either x or n must be specified.")
  if (!is.null(x)) n <- nrow(x)
  rf1 <- rep(NA, nf)
  rf2 <- rf1
  af <- rf1
  if (is.null(free1)) {
    #free1 <- simstandard::fixed2free(pop.mod1)
    pt <- lavaan::lavaanify(pop.mod1)
    pt <- pt[which(pt$op == "=~"),]
    free1 <- vector("list", length(unique(pt$lhs)))
    names(free1) <- unique(pt$lhs)
    for (i in unique(pt$lhs)) {
      free1[[i]] <- paste0(i, " =~ ", paste0(pt[which(pt$lhs == i), "rhs"], collapse = " + "))
    }
    free1 <- paste0(unlist(free1), collapse = "\n")
  }
  y <- try(lavaan::simulateData(
    model = pop.mod1,
    model.type = "cfa",
    sample.nobs = n,
    skewness = s,
    kurtosis = k,
    seed = seed,
    auto.fix.first = FALSE,
    std.lv = TRUE
  ),
  silent = TRUE)
  if (is.data.frame(y) & mode != "single") {
    sf1 <- try(lavaan::cfa(
      free1,
      y,
      estimator = "MLM",
      auto.fix.first = FALSE,
      std.lv = TRUE,
      warn = FALSE
    ),
    silent = TRUE)
    sf2 <- try(lavaan::cfa(
      free2,
      y,
      estimator = "MLM",
      auto.fix.first = FALSE,
      std.lv = TRUE,
      warn = FALSE
    ),
    silent = TRUE)
    if (inherits(sf1, "try-error") & inherits(sf2, "try-error")) {
      sf1 <- NA
      sf2 <- NA
    }
    if (typeof(sf1) == "S4" & typeof(sf2) == "S4") {
      if (sf1@Fit@converged & sf2@Fit@converged) {
        rf1 <- lavaan::fitmeasures(sf1)
        rf2 <- lavaan::fitmeasures(sf2)
      }
    }
  }
  if (is.data.frame(y) & mode == "single") {
    sf <- try(lavaan::cfa(
      free1,
      y,
      estimator = "MLM",
      auto.fix.first = FALSE,
      std.lv = TRUE,
      warn = FALSE
    ),
    silent = TRUE)
    if (inherits(sf, "try-error")) {
      sf <- NA
    }
    if (typeof(sf) == "S4") {
      if (sf@Fit@converged) {
        af <- lavaan::fitmeasures(sf)
      }
    }
  }
  if (mode != "single") {
    r <- rbind(rf1, rf2)
  }
  if (mode == "single") {
    r <- af
  }
  return(r)
}
