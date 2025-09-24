#' Helper function for structural model misspecification
#' @param pop.mod see gen_fit2
#' @param mc number of misspecified correlations
#' @param tc magnitude of misspecified correlations
#' @param seed ...
#' @keywords internal
#' @noRd
sm.misspec <- function(pop.mod,
                       mc = 1,
                       tc = NULL,
                       seed = NULL) {
  #mc: How many factor correlations should be misspecified?
  if (mc > 3)
    stop("More than three omitted correlations. So far, this is not possible.")
  #tc: Target correlation for misspecification, defaults to 0
  #Error handling
  if (is.null(tc))
    tc <- 0
  if (tc >= 1 |
      tc <= -1)
    stop("Please provide a proper correlation larger than -1 and smaller than 1")
  #Remove unintended whitespaces first for compatibility
  pop.mod <- stringr::str_replace_all(pop.mod, stringr::fixed(" "), "")
  #Make pop.mod a partable
  pt <- lavaan::lavaanify(pop.mod)
  #What are the factors?
  n.f <-
    unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(lhs))$lhs
  #Randomly pick factor correlations
  cs <- combn(n.f, 2)
  if (!is.null(seed))
    set.seed(seed)
  r.c <- matrix(cs[, sample(1:ncol(cs), mc, replace = F)], nrow = 2)
  #Find that factor correlation in pop.mod
  pmo <- strsplit(pop.mod, split = "\n")[[1]]
  sm <- rep(NA, ncol(r.c))
  #This is a bit tricky since both factors can be on each side of the ~~
  for (i in 1:ncol(r.c)) {
    in1 <-
      intersect(which(stringr::str_detect(pmo, pattern = paste0(
        as.character(r.c[1, i]), "~~"
      ))),
      which(stringr::str_detect(pmo, pattern = paste0(
        "\\*", as.character(r.c[2, i])
      ))))
    in2 <-
      intersect(which(stringr::str_detect(pmo, pattern = paste0(
        as.character(r.c[2, i]), "~~"
      ))),
      which(stringr::str_detect(pmo, pattern = paste0(
        "\\*", as.character(r.c[1, i])
      ))))
    #Replace the correlation by the new one
    if (length(in1) == 1) {
      sm[i] <- in1
      pmo[sm[i]] <- paste0(r.c[1, i], "~~", tc, "*", r.c[2, i])
    }
    if (length(in2) == 1) {
      sm[i] <- in2
      pmo[sm[i]] <- paste0(r.c[2, i], "~~", tc, "*", r.c[1, i])
    }
  }
  pop.mod.new <- paste0(pmo, collapse = "\n")
  return(pop.mod.new)
}
