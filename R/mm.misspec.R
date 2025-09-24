#' Helper function to misspecify population model with cross-loadings
#' @keywords internal
#' @param pop.mod see gen_fit2
#' @param ml number of cross-loadings
#' @param magni magnitude of cross-loadings
#' @param seed ...
#' @noRd
mm.misspec <- function(pop.mod,
                       ml = 1,
                       magni = NULL,
                       seed = NULL) {
  #Error handling
  if (ml > 3)
    stop("More than three cross-loadings specified. So far, this is not possible.")
  if (!is.null(magni)) {
    if (magni >= .9 |
        magni < 0.3)
      stop(
        "Please provide a proper (standardized) cross-loading magnitude below .9 and above .3."
      )
  }
  #ml: Number of misspecified loadings
  #Remove unintended whitespaces first for compatibility
  pop.mod <- stringr::str_replace_all(pop.mod, stringr::fixed(" "), "")
  #Make pop.mod a partable
  pt <- lavaan::lavaanify(pop.mod)
  #What are the factors?
  n.f <-
    unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(lhs))$lhs
  #What are the observed variables?
  n.ov <-
    unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(rhs))$rhs
  #Randomly pick an observed variable
  if (!is.null(seed))
    set.seed(seed)
  r.ov <- paste0("x", sample(1:length(n.ov), ml, replace = F))
  #Determine magnitude if not available for later
  if (is.null(magni)) {
    magni <-
      median((pt %>% dplyr::filter(op == "=~" &
                                     rhs %in% n.ov))$ustart, na.rm = T)
  }
  #What factor do they belong to?
  r.f <- rep(NA, ml)
  r.f <- vector("list", ml)
  for (i in 1:ml) {
    r.f[[i]] <-
      (pt %>% dplyr::filter(op == "=~" &
                              rhs %in% r.ov[i]) %>% dplyr::select(lhs))$lhs
    #Only first factor is taken when there are cross-loadings already (e.g., HB)
  }
  #Randomly pick a factor where r.ov is not a part of
  r.nf <- rep(NA, ml)
  for (i in 1:ml) {
    nn.f <- n.f[-which(n.f %in% r.f[[i]])]
    if (!is.null(seed))
      set.seed(seed)
    r.nf[i] <- sample(nn.f, 1, replace = F)
  }
  r.ov <- as_tibble(cbind(r.ov, r.nf))
  #Write the new meas. models
  mm <- vector("list", length(n.f))
  for (i in 1:length(n.f)) {
    lam <-
      (pt %>% dplyr::filter(op == "=~" &
                              rhs %in% n.ov &
                              lhs == n.f[i]))$ustart
    vs <-
      (pt %>% dplyr::filter(lhs == n.f[i] &
                              op == "=~") %>% dplyr::select(rhs))$rhs
    if (n.f[i] %in% r.nf) {
      lam <- c(lam, magni)
      vs <- c(vs, (r.ov %>% dplyr::filter(r.nf == n.f[i]))$r.ov)
    }
    tz <- rep(NA, length(lam))
    for (j in 1:length(lam)) {
      tz[j] <- paste0(lam[j], "*", vs[j])
    }
    mm[[i]] <- paste0(n.f[i], "=~", paste0(tz, collapse = "+"))
  }
  #Replace the pop.mod meas. models with the new ones
  pmo <- strsplit(pop.mod, split = "\n")[[1]]
  for (i in 1:length(n.f)) {
    pmo[which(stringr::str_detect(pmo, pattern = paste0(n.f[i], "=~")))] <-
      mm[[i]]
  }
  #Make this a pop.mod as before
  pop.mod.new <- paste0(pmo, collapse = "\n")
  return(pop.mod.new)
}
