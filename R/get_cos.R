#' Helpers to get cutoffs for different decision rules
#' @param co cutoff
#' @param tab table of quantiles
#' @param ind fit index
#' @param correct.fits see flex_co2
#' @param miss.fits ...
#' @keywords internal
#' @noRd
get_FCO1 <- function(co, tab, ind) {
fco1 <-
  as.numeric(
    tab %>% dplyr::filter(mod == "correct" & quant == co$alpha) %>%
      dplyr::select(dplyr::contains(ind))
  )
return(fco1)
}

get_FCO2 <- function(co, tab, ind) {
  ix <- index_guess(ind)
  if (ix == "GoF") {
    fco2 <- as.numeric(
      dplyr::if_else(
        tab %>% dplyr::filter(mod == "correct" & quant == co$alpha) %>%
          dplyr::select(dplyr::contains(ind)) > tab %>% dplyr::filter(mod == "miss" &
                                                                 quant == co$beta) %>%
          dplyr::select(dplyr::contains(ind)),
        tab %>% dplyr::filter(mod == "miss" &
                                quant == co$beta) %>%
          dplyr::select(dplyr::contains(ind)),
        tab %>% dplyr::filter(mod == "correct" &
                                quant == co$alpha) %>%
          dplyr::select(dplyr::contains(ind))
      )
    )
  }
  if (ix == "BoF") {
    fco2 <- as.numeric(
      dplyr::if_else(
        tab %>% dplyr::filter(mod == "correct" & quant == co$alpha) %>%
          dplyr::select(dplyr::contains(ind)) < tab %>% dplyr::filter(mod == "miss" &
                                                                 quant == co$beta) %>%
          dplyr::select(dplyr::contains(ind)),
        tab %>% dplyr::filter(mod == "miss" &
                                quant == co$beta) %>%
          dplyr::select(dplyr::contains(ind)),
        tab %>% dplyr::filter(mod == "correct" &
                                quant == co$alpha) %>%
          dplyr::select(dplyr::contains(ind))
      )
    )
  }
  return(fco2)
}

get_DFI <- function(co, tab, ind) {
  ix <- index_guess(ind)
  if (ix == "GoF") {
    dfi <- as.numeric(
      dplyr::if_else(
        tab %>% dplyr::filter(mod == "correct" & quant == co$alpha) %>%
          dplyr::select(dplyr::contains(ind)) > tab %>% dplyr::filter(mod == "miss" &
                                                                 quant == co$beta) %>%
          dplyr::select(dplyr::contains(ind)),
        tab %>% dplyr::filter(mod == "miss" &
                                quant == co$beta) %>%
          dplyr::select(dplyr::contains(ind)),
        NA
      )
    )
  }
  if (ix == "BoF") {
    dfi <- as.numeric(
      dplyr::if_else(
        tab %>% dplyr::filter(mod == "correct" & quant == co$alpha) %>%
          dplyr::select(dplyr::contains(ind)) < tab %>% dplyr::filter(mod == "miss" &
                                                                 quant == co$beta) %>%
          dplyr::select(dplyr::contains(ind)),
        tab %>% dplyr::filter(mod == "miss" &
                                quant == co$beta) %>%
          dplyr::select(dplyr::contains(ind)),
        NA
      )
    )
  }
  return(dfi)
}

get_Fix <- function(co, tab, ind) {
  #Single cutoffs only
  ind <- toupper(ind)
  if (ind == "SRMR")
    fix <- .08
  if (ind == "CFI")
    fix <- .95
  if (ind == "RMSEA")
    fix <- .06
  return(fix)
}

get_CP <- function(correct.fits, miss.fits, ind) {
  #Cutoff points from cutpointr
  ix <- index_guess(ind)
  dq <-
    c(dplyr::pull(correct.fits[, tolower(ind)]), dplyr::pull(miss.fits[, tolower(ind)]))
  fac <-
    c(rep("correct", nrow(correct.fits)), rep("miss", nrow(miss.fits)))
  dq <- dplyr::as_tibble(data.frame(ind = dq, fac = fac))
  cp <-  cutpointr::cutpointr(
    data = dq,
    x = ind,
    class = fac,
    direction = dplyr::if_else(ix == "GoF", ">=", "<="),
    na.rm = TRUE,
    silent = TRUE,
    method = cutpointr::maximize_metric,
    metric = cutpointr::F1_score
  )
  return(cp)
}

get_co <- function(co,
                   tab,
                   ind,
                   correct.fits = NULL,
                   miss.fits = NULL) {
  ct <- rep(NA, nrow(co))
  for (i in 1:nrow(co)) {
    if (co[i, "apr"] == "FCO1")
      ct[i] <- get_FCO1(co[i, ], tab, ind)
    if (co[i, "apr"] == "FCO2")
      ct[i] <- get_FCO2(co[i, ], tab, ind)
    if (co[i, "apr"] == "DFI")
      ct[i] <- get_DFI(co[i, ], tab, ind)
    if (co[i, "apr"] == "Fix")
      ct[i] <- get_Fix(co[i, ], tab, ind)
    if (co[i, "apr"] == "CP")
      ct[i] <- NA
  }
  res <- unname(ct)
  return(res)
}
