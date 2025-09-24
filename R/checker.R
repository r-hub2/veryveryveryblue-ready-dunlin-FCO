#' Helper function to check multiple issues with model and data.
#' @param pop.mod see gen_fit2
#' @param cfa ...
#' @param x ...
#' @param esti ...
#' @param sm ...
#' @param mm ...
#' @param rcv ...
#' @keywords internal
#' @noRd
checker <-
  function(pop.mod = pop.mod,
           cfa = cfa,
           x = x,
           esti = esti,
           sm = sm,
           mm = mm,
           rcv = rcv,
           cores = cores) {
    #Check cores
    if (parallel::detectCores() < cores) {
      warning(paste0("Fewer cores (", parallel::detectCores(), ") detected than specified (", cores, "). Please consider revision."))
    }
    #Is x a data.frame?
    if (!is.data.frame(x))
      stop("Argument x is not a data.frame. Please revise.")
    #Is the sample size larger than the number of moments (elements of the covariance matrix)?
    mod <- get_free_mod(pop.mod)
    pt <- lavaan::lavaanify(mod)
    x <- dplyr::as_tibble(x)
    n.ov <-
      unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(rhs))$rhs
    if (nrow(x) <= (length(n.ov) * (length(n.ov) - 1) / 2)) {
      warning("The dataset is smaller than the number of moments.")
    }
    #Is the manipulation permissible? More than one factor for sm / mm?
    n.f <-
      unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(lhs))$lhs
    if (length(n.f) == 1 &
        sm > 0 |
        length(n.f) == 1 &
        mm > 0) {
      stop("SM and MM manipulation cannot be done when the number of factors is 1.")
    }
    #Is sm & mm & rcv 0?
    if (sm == 0 & mm == 0 & rcv == 0)
      stop("No '000 model' can be used as the misspecified model.")
    #Is sm, mm, rcv max 3?
    if (sm > 3)
      stop("So far, no more than three misspecified correlations are supported.")
    if (mm > 3)
      stop("So far, no more than three misspecified loadings are supported.")
    if (rcv > 3)
      stop("So far, no more than three residual covariances are supported.")
    #Can there be so many restricted correlations?
    if (sm > (length(n.f) * (length(n.f) - 1) / 2)) {
      stop("More correlations are restricted than the model has.")
    }
    #Can there be so many cross-loadings?
    if (mm > length(n.f)) {
      stop("More cross-loadings are restricted than the model has. ")
    }
    #Can there be so many residual covariances?
    if (rcv > (length(n.ov) * (length(n.ov) - 1) / 2)) {
      stop("More residual covariances are considered than the model allows for.")
    }
    #Is the data mv-normal?
    sk <-
      median(psych::skew(x[, which(names(x) %in% n.ov)], na.rm = T), na.rm = T)
    ku <-
      median(psych::kurtosi(x[, which(names(x) %in% n.ov)], na.rm = T), na.rm = T)
    if (abs(sk) > 3 |
        abs(ku) > 7) {
      warning(
        paste0(
          "Skewness (",
          round(sk, 1),
          ") or Kurtosis (",
          round(ku, 1),
          ") are very high. Consider changing estimator and/or arguments sk and ku."
        )
      )
    }
    #Is the empirical covariance matrix positive definite?
    if (any(eigen(cov(x[, which(names(x) %in% n.ov)]))$values <= 0)) {
      stop("The empirical covariance matrix seems not to be positive definite.")
    }
    #Can the model be estimated?
    if (cfa) {
      test.res <-
        try(lavaan::cfa(model = mod,
                        data = x,
                        estimator = esti), silent = T)
    }
    if (!cfa) {
      test.res <-
        try(lavaan::sem(model = mod,
                        data = x,
                        estimator = esti), silent = T)
    }
    if (inherits(test.res, "try-error")) {
      stop("The model cannot be estimated, please revise.")
    }
    #Are dfs > 0?
    if (lavaan::lavInspect(test.res, "fit")["df"] <= 0) {
      stop("The implied model has no positive degrees of freedom.")
    }
    #Did the model converge?
    if (!lavaan::lavInspect(test.res, "converged")) {
      stop("The model has not converged.")
    }
    #Are there negative variances or non-positive definite matrices?
    if (any(eigen(lavaan::lavInspect(test.res, "cov.ov"))$values <= 0)) {
      stop("The implied covariance matrix of the indicators seems not be positive definite.")
    }
    if (any(eigen(lavaan::lavInspect(test.res, "cov.lv"))$values <= 0)) {
      stop(
        "The implied covariance matrix of the latent variables seems not be positive definite."
      )
    }
    if (any(diag(lavaan::lavInspect(test.res, "theta")) <= 0)) {
      warning("The model has implied some negative variances.")
    }
  }

#Combine
build_and_check <- function(seed, sm, mm, rcv, pop.mod, cfa, x, esti, cores) {
  mpop <- model.builder(seed = seed, sm = sm, mm = mm, rcv = rcv, pop.mod = pop.mod)
  checker(pop.mod = mpop, cfa = cfa, x = x, esti = esti, sm = sm, mm = mm, rcv = rcv, cores = cores)
  mpop
}
