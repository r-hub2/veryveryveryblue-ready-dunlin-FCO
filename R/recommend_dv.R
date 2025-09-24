#' Obtain recommendations for discriminant validity testing
#'
#' This function recommends on potential issues for discriminant validity testing, based on  differences between fit values and differences between flexible cutoffs. Two approaches of testing are supported: merging and constraining.
#' @param fits A list of simulated fit indices obtained from gen_fit. Based on the structure of fits, the number of models is derived.
#' @param index A vector of fit indices or measures provided by function fitmeasures in package lavaan. The default is set to CFI.
#' @param digits An optional integer to round fit values and cutoffs (min: 1, max: 5).
#' @return A list of information regarding discriminant validity testing.
#' @examples
#' \donttest{
#'#Note: Demonstration only! Please use higher numbers of replications for your applications (>= 500).
#'mod <- "
#'F1 =~ Q5 + Q7 + Q8
#'F2 =~ Q2 + Q4
#'F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
#'F4 =~ Q1 + Q17
#'F5 =~ Q6 + Q14 + Q15 + Q16
#'"

#'#Two models for discriminant validity testing, this resembles constraining with a cutoff of .9
#'fits.dv.con <- gen_fit(
#'  mod1 = mod,
#'  x = bb1992,
#'  rep = 10,
#'  dv = TRUE,
#'  dv.factors = c("F4", "F5"),
#'  dv.cutoff = .9
#')
#'recommend_dv(fits.dv.con)
#' #Two models for discriminant validity testing, this resembles merging.
#'fits.dv.merge <- gen_fit(
#'  mod1 = mod,
#'  x = bb1992,
#'  rep = 10,
#'  dv = TRUE,
#'  dv.factors = c("F4", "F5"),
#'  merge.mod = TRUE
#')
#'recommend_dv(fits.dv.merge)
#'}
#' @export
recommend_dv <-
  function(fits,
           index = "CFI",
           digits = 3) {
    if (is.null(dim(fits$fco[[1]])))
      stop(
        "Only dual model fit indices are supported so far for this function. Please revise or use flex_co."
      )
    checkmate::assertCharacter(fits$mod1,
                               fixed = "=~")
    checkmate::assertDataFrame(
      fits$x,
      min.rows = 50,
      min.cols = 4,
      col.names = "unique"
    )
    checkmate::assertVector(fits$dv.factors, len = 2, null.ok = TRUE)
    checkmate::assertLogical(fits$merge.mod, null.ok = TRUE)
    checkmate::assertVector(fits$dv.cutoff, len = 1, null.ok = TRUE)
    checkmate::assert_int(digits, lower = 1, upper = 5)
    mode <- ifelse(fits$merge.mod, "merging", "constraining")
    fm <-
      try(lavaan::fitmeasures(
        lavaan::cfa(
          fits$mod1,
          data = fits$x,
          estimator = "MLM",
          auto.fix.first = FALSE,
          std.lv = TRUE
        )
      ), silent = TRUE)
    if (inherits(fm, "try-error"))
      stop("Invalid model or data. Please revise.")
    # if (!tolower(index) %in% names(fm))
    #   stop(
    #     "The index names provided do not match the names of the indices provided by lavaan. Please revise."
    #   )
    if (mode == "constraining") {
      mod2 <-
        constr_mod(fits$mod1,
                   dv.factors = fits$dv.factors,
                   dv.cutoff = fits$dv.cutoff)
      fm2 <-
        try(lavaan::fitmeasures(
          lavaan::cfa(
            fits$mod2,
            data = fits$x,
            estimator = "MLM",
            auto.fix.first = FALSE,
            std.lv = TRUE
          )
        ), silent = TRUE)
      if (inherits(fm2, "try-error"))
        stop("Invalid model or data. Please revise.")
    }
    if (mode == "merging") {
      pop.mod <- pop_mod(
        mod = fits$mod1,
        x = fits$x,
        type = fits$type,
        standardized = fits$standardized
      )$pop.mod
      pt2 <-
        try(merge_factors(lavaan::cfa(pop.mod, fits$x, warn = FALSE),
                          merged.factors = fits$dv.factors),
            silent = TRUE)
      if (inherits(pt2, "try-error"))
        stop("Merging two-factors not successful. Please check.")
      mod2 <- get_free(pt2, fits$dv.factors, mode)
      fm2 <-
        try(lavaan::fitmeasures(
          lavaan::cfa(
            mod2,
            data = fits$x,
            estimator = "MLM",
            auto.fix.first = FALSE,
            std.lv = TRUE
          )
        ), silent = TRUE)
      if (inherits(fm2, "try-error"))
        stop("Invalid model or data. Please revise.")
    }
    index <- tolower(index)
    ap <- c(.001, .01, .05, .10)
    #gof <- ifelse(sapply(index, index_guess) == "GoF", TRUE, FALSE)
    fc <-
      suppressWarnings(sapply(ap, flex_co, fits = fits, index = index))
    #fc <- suppressWarnings(sapply(ap, flex_co, fits = fits, index = index, gof = gof))
    #fc2 <- suppressWarnings(sapply(ap, flex_co, fits = fits, index = index, gof = !gof))
    #for (i in seq_len(length(ap))) {
    #  fc[1,][[i]][2, ] <- fc2[1,][[i]][2, ]
    #}
    co <- lapply(fc[1,], as.data.frame)
    gof <- fc[5,]
    na <- fc[7, 1]
    sh.na <- fc[8, 1]
    rf1 <- fm[tolower(index)]
    rf2 <- fm2[tolower(index)]
    tab <- data.table::rbindlist(co)
    mn <- ifelse(mode == "constraining", "constrained", "merged")
    tab <-
      as.data.frame(cbind(tab, expand.grid(c("original", mn), ap)))
    colnames(tab) <- c(toupper(index), "model", "alpha")
    fi <- data.frame(rbind(rf1, rf2))
    rownames(fi) <- NULL
    alpha <- NULL
    fi$model <- c("original", mn)
    colnames(fi)[seq_len(length(index))] <- c(toupper(index))
    diff <-
      data.frame(matrix(NA, nrow = length(ap), ncol = length(index)))
    for (i in seq_len(length(ap))) {
      ss <- subset(tab, alpha == ap[i])[, seq_len(length(index))]
      if (length(index) == 1) {
        diff[i,] <- ss[1] - ss[2]
      }
      if (length(index) > 1) {
        diff[i, ] <- ss[1, ] - ss[2, ]
      }
    }
    names(diff) <- toupper(index)
    diff <-
      t(rbind(diff, fi[1, seq_len(length(index))] - fi[2, seq_len(length(index))]))
    colnames(diff) <- c(paste0("cutoff ", ap), "fit")
    decide.dv <- function(tab, fi, diff, index) {
      model <- NULL
      ctab <- subset(tab, model == "original")[, toupper(index)]
      cfi <- subset(fi, model == "original")[, toupper(index)]
      #whr <- which(sapply(colnames(ctab), index_guess) == "BoF")
      whr <- which(vapply(colnames(ctab), index_guess, character(1)) == "BoF")
      if (!is.null(nrow(ctab))) {
        decs <- data.frame(matrix(NA, nrow = nrow(ctab), ncol = ncol(ctab)))
        for (i in seq_len(nrow(decs))) {
          #For GoF, that is the version from the flc paper
          decs[i, ] <- ifelse(cfi <= ctab[i, ], 1,-1)
        }
        #For BoF: reverse
        decs[, whr] <- decs[, whr] * -1
        decs[decs == 1] <- "confirmed"
        decs[decs == -1] <- "rejected"
      }
      if (is.null(nrow(ctab))) {
        #For GoF, that is the version from the flc paper
        decs <- ifelse(cfi <= ctab, 1,-1)
        #For BoF: reverse
        decs[whr] <- decs[whr] * -1
        decs[decs == 1] <- "confirmed"
        decs[decs == -1] <- "rejected"
      }
      return(decs)
    }
    # decide.dv <- function(diff) {
    #   whc <- which(colnames(diff) == "fit")
    #   whr <- which(sapply(rownames(diff), index_guess) == "BoF")
    #   whr <- which(vapply(colnames(ctab), index_guess, character(1)) == "BoF")
    #   decs <- data.frame(matrix(NA, nrow = nrow(diff), ncol = ncol(diff) - 1))
    #   for (i in seq_len(ncol(decs))) {
    #     #For GoF, that is the version from the flc paper
    #     decs[,i] <- ifelse(diff[,whc] >= diff[,i], 1, -1)
    #   }
    #   #For BoF: reverse
    #   decs[whr,] <- decs[whr,] * -1
    #   decs[decs == 1] <- "confirmed"
    #   decs[decs == -1] <- "rejected"
    #   return(decs)
    # }
    #decs <- decide.dv(diff)
    decs <- decide.dv(tab, fi, diff, index)
    #decs <- data.frame(t(decs))
    diff <- data.frame(t(diff))
    decs <- data.frame(decs)
    names(decs) <- colnames(diff)
    rownames(decs) <- rownames(diff)[seq_len(length(ap))]
    tab[, toupper(index)] <- round(tab[, toupper(index)], digits)
    fi[, toupper(index)] <- round(fi[, toupper(index)], digits)
    diff[, toupper(index)] <- round(diff[, toupper(index)], digits)
    if (fits$rep < 500)
      warning(
        "The number of replications is lower than the recommended minimum of 500. Consider with care."
      )
    res <- list(
      "cutoffs" = tab,
      "fit.values" = fi,
      "differences" = diff,
      "decisions" = decs,
      "replications" = fits$rep,
      "comment" = paste0(
        "Approach for discriminant validity testing: ",
        mode,
        ". Discriminant validity is confirmed if the fit value from the constrained/merged model is smaller (GoF) / larger (BoF) than the respective cutoff of original model."
      )
    )
    return(res)
  }
