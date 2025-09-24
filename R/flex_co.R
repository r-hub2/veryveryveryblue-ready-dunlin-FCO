#' Obtain flexible cutoffs for one or two models
#'
#' @param fits A list of simulated fit indices obtained from gen_fit. Based on the structure of fits, the number of models is derived.
#' @param index A vector of fit indices or measures provided by function fitmeasures in package lavaan
#' @param alpha.lev The predefined uncertainty
#' For example, if the default uncertainty of .05 (5 percent) is accepted a-priori, the 5 percent stats::quantile (of type 8, see ?stats::quantile) of the simulated distribution for correctly specified CFA models with the given model and sample characteristics determines the flexible cutoff.
#' Options are .001, .01, .05, and .10. Higher values are more conservative.
#' @param gof An optional vector as to whether the indices are GoF (Goodness-of-Fit index)? If TRUE, a GoF is assumed. If FALSE, a BoF is assumed.
#' Depending on the nature of the underlying fit index, the appropriate lower (GoF) or upper (BoF) width of the respective confidence interval as defined by the stats::quantile is used to derive the flexible cutoff.
#' If not provided or not equal to the number of fit indices, the function guesses the type for known fit indices (e.g., SRMR is a BoF).
#' @return A list of information regarding the selected fit index providing its flexible cutoff for the given parameters.
#' @examples
#'#Note: Demonstration only! Please use higher numbers of replications for your applications (>= 500).
#'#A single model to obtain fit indices for
#'mod <- "
#'F1 =~ Q5 + Q7 + Q8
#'F2 =~ Q2 + Q4
#'F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
#'F4 =~ Q1 + Q17
#'F5 =~ Q6 + Q14 + Q15 + Q16
#'"
#'fits.single <- gen_fit(mod1 = mod, x = bb1992, rep = 10, standardized = FALSE)
#'flex_co(fits = fits.single, index = c("CFI", "SRMR"))
#'\donttest{
#' #Two models, an unconstrained and a constrained model to compare fit indices
#'mod.con <- "
#'F1 =~ Q5 + Q7 + Q8
#'F2 =~ Q2 + Q4
#'F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
#'F4 =~ Q1 + Q17
#'F5 =~ Q6 + Q14 + Q15 + Q16
#'F1 ~~ 0 * F2
#'"
#'fits.con <- gen_fit(
#'  mod1 = mod,
#'  mod2 = mod.con,
#'  x = bb1992,
#'  rep = 10
#')
#'flex_co(fits = fits.con,
#'        index = c("CFI", "SRMR"),
#'        alpha.lev = .05)
#'
#' #Two models for discriminant validity testing, this resembles constraining with a cutoff of .9
#'fits.dv.con <- gen_fit(
#'  mod1 = mod,
#'  x = bb1992,
#'  rep = 10,
#'  dv = TRUE,
#'  dv.factors = c("F4", "F5"),
#'  dv.cutoff = .9
#')
#'flex_co(fits = fits.dv.con,
#'index = "CFI",
#'alpha.lev = .05)
#'mod.dv.con <- "
#'F1 =~ Q5 + Q7 + Q8
#'F2 =~ Q2 + Q4
#'F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
#'F4 =~ Q1 + Q17
#'F5 =~ Q6 + Q14 + Q15 + Q16
#'F4 ~~ .9 * F5
#'"
#'lavaan::fitmeasures(
#'  lavaan::cfa(
#'    model = mod.dv.con,
#'    data = bb1992,
#'    auto.fix.first = FALSE,
#'    std.lv = TRUE
#'  ),
#'  fit.measures = "cfi"
#')

#' #Two models for discriminant validity testing, this resembles merging.
#'fits.dv.merge <- gen_fit(
#'  mod1 = mod,
#'  x = bb1992,
#'  rep = 10,
#'  dv = TRUE,
#'  dv.factors = c("F4", "F5"),
#'  merge.mod = TRUE)
#'
#'flex_co(fits = fits.dv.merge,
#'index = "CFI",
#'alpha.lev = .05)
#'mod.dv.merge <- "
#'F1 =~ Q5 + Q7 + Q8
#'F2 =~ Q2 + Q4
#'F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
#'F4 =~ Q1 + Q17 + Q6 + Q14 + Q15 + Q16
#'"
#'lavaan::fitmeasures(
#'  lavaan::cfa(
#'    model = mod.dv.merge,
#'    data = bb1992
#'  ),
#'  fit.measures = "cfi"
#')
#'}
#' @export
flex_co <-
  function(fits,
           index,
           alpha.lev = .05,
           gof = NULL) {
    checkmate::assertCharacter(index, min.chars = 3)
    index <- tolower(index)
    checkmate::assertNumeric(alpha.lev, len = 1)
    checkmate::assertLogical(gof, null.ok = TRUE)
    if (length(fits$fco) < 500)
      warning(
        "The number of replications is lower than the recommended minimum of 500. Consider with care."
      )
    if (!alpha.lev %in% c(.001, .01, .05, .10))
      warning("The selected alpha level is unequal to .001, .01, .05, or .10.")
    if (is.null(gof) | length(gof) != length(index)) {
      gof <- ifelse(vapply(index, index_guess, character(1)) == "GoF", TRUE, FALSE)
    }
    ml <- ifelse(is.null(nrow(fits$fco[[1]])), 1, 2)
    fits <- fits$fco
    if (ml == 1) {
      #fits.nna <- fits[!sapply(fits, function(x)
      #  all(is.na(x)))]
      fits.nna <- fits[!vapply(fits, function(x)
        all(is.na(x)), numeric(1))]
      na <- length(fits) - length(fits.nna)
      sh.na <- na / length(fits)
      if (!all(index %in% names(fits.nna[[1]])))
        stop("At least one selected index is not a supported fitmeasure in lavaan.")
      #sf <- unname(sapply(fits.nna, function(x)
      #  x[index]))
      sf <- unname(vapply(fits.nna, function(x)
        x[index], numeric(length(index))))
      co <- rep(NA, length(index))
      #probs <-
      #  sapply(gof, function(x)
      #    ifelse(x, alpha.lev, 1 - alpha.lev))
      probs <-
        vapply(gof, function(x)
          ifelse(x, alpha.lev, 1 - alpha.lev), numeric(1))
      for (i in seq_len(length(index))) {
        #Type of stats::quantile is set to 8
        if (length(index) > 1) {
          co[i] <-
            stats::quantile(sf[i,], probs[i], type = 8)
        }
        if (length(index) == 1) {
          co[i] <-
            stats::quantile(sf, probs, type = 8)
        }
      }
      names(co) <- toupper(index)
      names(gof) <- toupper(index)
      rco <-
        list(
          "cutoff" = co,
          "index" = toupper(index),
          "alpha" = alpha.lev,
          "gof" = gof,
          "replications" = length(fits),
          "number of non-converging models" = na,
          "share of non-converging models" = round(sh.na, 3)
        )
    }
    if (ml == 2) {
      #fits.nna <- fits[!sapply(fits, function(x)
      #  all(is.na(x)))]
      fits.nna <- fits[!vapply(fits, function(x)
        all(is.na(x)), numeric(1))]
      na <- length(fits) - length(fits.nna)
      sh.na <- na / length(fits)
      if (!all(index %in% names(fits.nna[[1]][1,])))
        stop("At least one selected index is not a supported fitmeasure in lavaan.")
      #vf1 <-
      #  unname(sapply(lapply(fits.nna, function(x)
      #    x[1,]), function(x)
      #      x[index]))
      vf1 <-
        unname(vapply(lapply(fits.nna, function(x)
          x[1,]), function(x)
            x[index], numeric(length(index))))
      #vf2 <-
      #  unname(sapply(lapply(fits.nna, function(x)
      #    x[2,]), function(x)
      #      x[index]))
      vf2 <-
        unname(vapply(lapply(fits.nna, function(x)
          x[2,]), function(x)
            x[index], numeric(length(index))))
      co1 <- rep(NA, length(index))
      co2 <- co1
      #probs <-
      #  sapply(gof, function(x)
      #    ifelse(x, alpha.lev, 1 - alpha.lev))
      probs <-
        vapply(gof, function(x)
          ifelse(x, alpha.lev, 1 - alpha.lev), numeric(1))
      for (i in seq_len(length(index))) {
        #Type of stats::quantile is set to 8
        if (length(index) > 1) {
          co1[i] <-
            stats::quantile(vf1[i,], probs[i], type = 8)
          co2[i] <-
            stats::quantile(vf2[i,], probs[i], type = 8)
        }
        if (length(index) == 1) {
          co1[i] <-
            stats::quantile(vf1, probs, type = 8)
          co2[i] <-
            stats::quantile(vf2, probs, type = 8)
        }
      }
      co <- rbind(co1, co2)
      dimnames(co) <- list(c("Model 1", "Model 2"), toupper(index))
      names(gof) <- toupper(index)
      rco <-
        list(
          "cutoff" = co,
          "difference" = co[1,] - co[2,],
          "index" = toupper(index),
          "alpha" = alpha.lev,
          "gof" = gof,
          "replications" = length(fits),
          "number of non-converging models" = na,
          "share of non-converging models" = round(sh.na, 3)
        )
    }
    return(rco)
  }
