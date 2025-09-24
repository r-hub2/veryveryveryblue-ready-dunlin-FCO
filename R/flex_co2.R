#' Obtain cutoffs for fit indices simulated by gen_fit2.
#' @details This function returns cutoffs for the the fit indices simulated by gen_fit2. For details, please refer to gen_fit2. Please note that the results are only based on the simulation of the model specified and its misspecified variant under the conditions of the gen_fit2 function.
#' @param fits The object returned from gen_fit2 which is a list of correct fit values and misspecified fit values.
#' @param correct.fits For compatibility reasons, the correct fit values can be defined separately.
#' @param miss.fits For compatibility reasons, the misspecified fit values can be defined separately.
#' @param index A vector of length >= 1 with names of fit indices the user wishes to explore. Capitalization does not matter, either e.g., CFI or cfi are accepted. Default is CFI as it might be easier to understand.
#' @param alpha The acceptable Type I error representing the empirical quantile p, see details in gen_fit2. Multiple values can be provided as a vector.  Default is .05.
#' @param beta The acceptable Type II error representing the empirical quantile p, see details in gen_fit2. Multiple values can be provided as a vector. Default is c(.05, .10).
#' @return A list consisting of a tibble for the empirical quantiles estimated for the provided alpha, beta values and indices, a tibble for the derived cutoff values for each parameter (alpha, beta, approach, index, cutoff), a tibble for the evaluation (also True Negatives, False Positives, True Positives, False Negatives, Type I error, Type II error, Sum of both Types, Power and Specificity), a vector for the notation on the evaluation tibble, and a tibble displaying overlap statistics for each index (Overlap percentage, AUC, U-test).
#' @examples
#' #Simple example
#' library(lavaan)
#' library(dplyr)
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'              speed   =~ x7 + x8 + x9 '
#'
#' fit <- cfa(
#'   HS.model,
#'   data = HolzingerSwineford1939
#' )
#' #Note: Demonstration only! Please use higher numbers of replications for your applications (>= 500).
#' fits <- gen_fit2(fit = fit, rep = 10)
#' #Default evaluation:
#' flex_co2(fits)
#' #Changed alpha and beta values:
#' flex_co2(fits, alpha = .05, beta = .05)
#' flex_co2(fits, alpha = .10, beta = .20)
#' #Different fit indices:
#' flex_co2(fits, index = c("CFI", "SRMR", "RMSEA"))
#' @export
flex_co2 <- function(fits = NULL,
         correct.fits = NULL,
         miss.fits = NULL,
         index = "CFI",
         alpha = c(.05),
         beta = c(.05, .10)) {
  if (is.null(fits) & is.null(correct.fits) & is.null(miss.fits)) {
    stop(
      "Please provide the list of fits returned by gen_fits2 or specify correct.fits and miss.fits accordingly."
    )
  }
  if (!is.null(fits)) {
    cf <- fits[[1]]
    mf <- fits[[2]]
  }
  if (is.null(fits)) {
    cf <- correct.fits
    mf <- miss.fits
  }
  ind <- toupper(index)
  #Determine quantiles
  cs <-
    as_tibble(matrix(NA, ncol = length(ind), nrow = length(alpha)))
  names(cs) <- ind
  for (i in 1:ncol(cs)) {
    if (index_guess(names(cs)[i]) == "GoF") {
      cs[, i] <- quantile(cf %>% dplyr::select(dplyr::any_of(tolower(ind[i]))),
                          alpha,
                          type = 8,
                          na.rm = TRUE)
    }
    if (index_guess(names(cs)[i]) == "BoF") {
      cs[, i] <- quantile(cf %>% dplyr::select(dplyr::any_of(tolower(ind[i]))),
                          1 - alpha,
                          type = 8,
                          na.rm = TRUE)
    }
  }
  ms <-
    as_tibble(matrix(NA, ncol = length(ind), nrow = length(beta)))
  names(ms) <- ind
  for (i in 1:ncol(ms)) {
    if (index_guess(names(ms)[i]) == "GoF") {
      ms[, i] <- quantile(mf %>% dplyr::select(dplyr::any_of(tolower(ind[i]))),
                          1 - beta,
                          type = 8,
                          na.rm = TRUE)
    }
    if (index_guess(names(ms)[i]) == "BoF") {
      ms[, i] <- quantile(mf %>% dplyr::select(dplyr::any_of(tolower(ind[i]))),
                          beta,
                          type = 8,
                          na.rm = TRUE)
    }
  }
  tab <- rbind(cs, ms)
  tab$quant <- c(alpha, beta)
  tab$mod <-
    c(rep("correct", length(alpha)), rep("miss", length(beta)))
  tabl <- tidyr::pivot_longer(tab, cols = 1:length(ind), names_to = "index")
  tabl$type <- sapply(tabl$index, index_guess)
  tabl <- tabl %>% dplyr::mutate(quant = dplyr::if_else(
    type == "GoF" &
      mod == "miss",
    1 - quant,
    dplyr::if_else(type == "BoF" &
      mod == "correct",
      1 - quant, quant),
    quant)) %>% dplyr::arrange(index, mod)
  #Cutoffs
  apr <- c("FCO1", "FCO2", "DFI", "Fix", "CP")
  co <- tidyr::expand_grid(alpha = alpha,
                           beta = beta,
                           apr = apr)
  co <-
    dplyr::bind_cols(co, matrix(
      NA,
      nrow = nrow(co),
      ncol = length(ind),
      dimnames = list(rep(NA, nrow(co)), ind)
    ))
  #Do get_CP manually and fill
  tmp <- vector("list", length(ind))
  names(tmp) <- ind
  for (i in ind) {
    co[, i] <- get_co(
      co,
      tab,
      ind = i,
      correct.fits = cf,
      miss.fits = mf
    )
    tmp[[i]] <- get_CP(correct.fits = cf,
                       miss.fits = mf,
                       ind = i)
    co[which(co$apr == "CP"), i] <- tmp[[which(names(tmp) == i)]]$optimal_cutpoint
  }
  #Power
  nop <-
    "Notation for power: positive = Models are different, negative = Models are equal, TN = correct & GoF > co or BoF < co, TP = miss & GoF < co or BoF > co, FP = correct & GoF < co or BoF > co, FN = miss & GoF > co or BoF < co"
  co <-
    co %>% tidyr::pivot_longer(cols = dplyr::contains(ind),
                               names_to = "index",
                               values_to = "cutoff")
  cp <- data.frame(matrix(NA, nrow = nrow(co), ncol = 9))
  names(cp) <-
    c("TN",
      "FP",
      "TP",
      "FN",
      "TypeI",
      "TypeII",
      "SumTypes",
      "Power",
      "Specificity")
  for (i in 1:nrow(co)) {
    cp[i, ] <-
      power.finder(cutoff = co$cutoff[i],
                   index = co$index[i],
                   cf,
                   mf)
  }
  cp <- dplyr::bind_cols(co, cp)
  cp <- cp %>% dplyr::arrange(SumTypes)
  #Overlap and tests
  ol <- dplyr::as_tibble(matrix(NA, nrow = 3, ncol = length(ind) + 1))
  names(ol) <- c("Statistic", ind)
  ol[, 1] <-
    c("Overlap (percentage)",
      "AUC (Area under curve)",
      "U-test (effect size d)")
  for (i in ind) {
    cft <-
      na.omit(unname(unlist(cf %>% dplyr::select(
        dplyr::all_of(tolower(i))
      ))))
    mft <-
      na.omit(unname(unlist(mf %>% dplyr::select(
        dplyr::all_of(tolower(i))
      ))))
    ol[1, i] <-
      overlapping::overlap(list(cft, mft), type = "2")$OV
    ol[2, i] <- tmp[[which(names(tmp) == i)]]$AUC
    wr <- rcompanion::wilcoxonR(x = c(cft, mft), g = factor(c(rep(
      "c", length(cft)
    ), rep(
      "m", length(mft)
    ))))
    ol[3, i] <- 2 * wr / sqrt(1 - wr^2)
  }
  co <- co %>%
    dplyr::rename(dec.rule = apr)
  cp <- cp %>%
    dplyr::rename(dec.rule = apr)
  res <-
    list(
      quantiles = tabl,
      cutoffs = co,
      evaluation = cp,
      notation = nop,
      overlap = ol
    )
  return(res)
}
