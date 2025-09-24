#' Helper function for residual covariance misspecification
#' @keywords internal
#' @param pop.mod see gen_fit2
#' @param me number of rcvs misspecified
#' @param te magnitude of misspecified rcvs
#' @param seed ...
#' @noRd
rv.misspec <- function(pop.mod,
                       me = 1,
                       te = NULL,
                       seed = NULL) {
  #me: Number of misspecified correlated residual variances
  #te: Target correlation of misspecified correlated residual variances, defaults to .3
  #Error handling
  if (me > 3)
    stop("More than three omitted residual correlations. So far, this is not possible.")
  if (is.null(te))
    te <- .3
  if (te >= 1 |
      te <= -1)
    stop("Please provide a proper residual correlation larger than -1 and smaller than 1")
  #Remove unintended whitespaces first for compatibility
  pop.mod <- stringr::str_replace_all(pop.mod, stringr::fixed(" "), "")
  #Make pop.mod a partable
  pt <- lavaan::lavaanify(pop.mod)
  #What are the observed variables?
  n.ov <-
    unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(rhs))$rhs
  #Randomly pick residual variance correlations
  #This variant makes sure that the correlations are from different pairs
  es <- matrix(NA, nrow = 2, ncol = me)
  for (i in 1:me) {
    if (i == 1) {
      if (!is.null(seed))
        set.seed(seed)
      es[, i] <- sample(n.ov, 2, replace = F)
      o.ov <- n.ov[!n.ov %in% es[, i]]
    }

    if (i > 1) {
      if (!is.null(seed))
        set.seed(seed)
      es[, i] <- sample(o.ov, 2, replace = F)
      o.ov <- o.ov[!o.ov %in% es[, i]]
    }
  }
  #Add residual variance correlations
  pmo <- strsplit(pop.mod, split = " \n")[[1]]
  pmo <- c(pmo, rep(NA, me))
  for (i in 1:me) {
    pmo[length(pmo) + i - me] <-
      paste0(es[1, i], "~~", te, "*", es[2, i])
  }
  pop.mod.new <- paste0(pmo, collapse = "\n")
  return(pop.mod.new)
}
