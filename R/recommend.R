#' Obtain recommendations based on Mai et al. (2021)
#'
#' This function recommends pre-defined selected fit indices in case the user does not know which fit index should be used for model evaluation. Results may differ based on three settings, the sample size of the data, the research purpose of the investigated model and the focus of the model. For obvious reasons, this function only works for single models and does not accept any other model type.
#' @param fits A list of simulated fit indices obtained from gen_fit. Based on the structure of fits, the number of models is derived.
#' @param purpose The research purpose of the model investigated. Is the underlying model novel (default) or established (= established). This parameter is relevant to find the proper recommended fit indices.
#' @param focus The focus of estimation for the model. Is the focus on CFA (default) or analyzing the structural model of a theoretical model (= structural)? This parameter is relevant to find the proper recommended fit indices.
#' @param override Should the recommendations by Mai et al. (2021) overridden (default: FALSE)?  This may be useful to explore models outside of the scope of the paper. In this case, the recommended fit indices are not determined by the function, and hence need to be provided.
#' In this case, the function requires the argument index.
#' @param index An optional vector of fit indices or measures provided by function fitmeasures in package lavaan. This argument is required when override is TRUE. It is ignored otherwise.
#' @param digits An optional integer to round fit values and cutoffs (min: 1, max: 5).
#' @return A list of information regarding the recommended fit indices based on Mai et al. (2021) or when overridden, based on the provided indices.
#' @references Mai, R., Niemand, T., & Kraus, S. (2021). A Tailor-Fit Model Evaluation Strategy for Better Decisions about Structural Equation Models, Technological Forecasting & Social Change, 173(December) 121142. https://doi.org/10.1016/j.techfore.2021.121142
#' @examples
#'#Note: Demonstration only! Please use higher numbers of replications for your applications (>= 500).
#'mod <- "
#'F1 =~ Q5 + Q7 + Q8
#'F2 =~ Q2 + Q4
#'F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
#'F4 =~ Q1 + Q17
#'F5 =~ Q6 + Q14 + Q15 + Q16
#'"
#'fits.single <- gen_fit(mod1 = mod, x = bb1992, rep = 10, standardized = FALSE)
#'recommend(fits.single)
#'recommend(fits.single, purpose = "established")
#'recommend(fits.single,
#'          override = TRUE,
#'          index = c("CFI", "SRMR"))
#' @export
recommend <-
  function(fits,
           purpose = "novel",
           focus = "cfa",
           override = FALSE,
           index = NULL,
           digits = 3) {
    if (!is.null(dim(fits$fco[[1]])))
      stop(
        "Only single model fit indices are supported so far for this function. Please revise or use flex_co."
      )
    if (is.null(fits$pop.mod1)) {
      checkmate::assertCharacter(fits$mod1,
                                 fixed = "=~")
      checkmate::assertDataFrame(
        fits$x,
        min.rows = 50,
        min.cols = 4,
        col.names = "unique"
      )
      n <- nrow(fits$x)
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
    }
    if (is.null(fits$mod1)) {
      checkmate::assertCharacter(fits$pop.mod1,
                                 fixed = "*")
      checkmate::assertNumeric(fits$n, lower = 50, upper = 50000, null.ok = TRUE)
      x <- lavaan::simulateData(model = fits$pop.mod1, sample.nobs = fits$n, seed = fits$seed)
      fm <-
        try(lavaan::fitmeasures(
          lavaan::cfa(
            simstandard::fixed2free(fits$pop.mod1),
            data = x,
            estimator = "MLM",
            auto.fix.first = FALSE,
            std.lv = TRUE
          )
        ), silent = TRUE)
      if (inherits(fm, "try-error"))
        stop("Invalid model or data. Please revise.")
    }
    if (purpose != "novel")
      purpose <- "established"
    checkmate::assert(
      checkmate::checkCharacter(purpose,
                                pattern = "novel"),
      checkmate::checkCharacter(purpose,
                                pattern = "established"),
    )
    focus <- tolower(focus)
    if (focus != "cfa")
      focus <- "structure"
    checkmate::assert(
      checkmate::checkCharacter(focus,
                                pattern = "cfa"),
      checkmate::checkCharacter(focus,
                                pattern = "structure"),
    )
    checkmate::assert_int(digits, lower = 1, upper = 5)
    if (!override) {
      if (purpose == "established" & focus == "cfa" & n <= 200) {
        #C1: SRMR flex
        index <- "srmr"
        ap <- c(.001, .01, .05, .10)
        #fc <-
        #  suppressWarnings(sapply(ap, flex_co, fits = fits, index = index))
        fc <- suppressWarnings(lapply(
          ap, flex_co, fits = fits, index = index
        ))
        co <- vector("list", length(fc))
        gof <- co
        for (i in seq_len(length(fc))) {
          co[[i]] <- fc[[i]]$cutoff
          gof[[i]] <- fc[[i]]$gof
        }
        na <- fc[[1]]$`number of non-converging models`
        sh.na <- fc[[1]]$`share of non-converging models`
        #co <- fc[1,]
        #gof <- fc[4,]
        #na <- fc[6, 1]
        #sh.na <- fc[7, 1]
        rf <- fm[index]
      }
      if (purpose == "established" & focus == "cfa" & n > 200) {
        #C2: CFI fix
        index <- "cfi"
        ap <- NULL
        fc <- NULL
        na <- NULL
        sh.na <- NULL
        co <- .95
        gof <- index_guess(index)
        rf <- fm[index]
      }
      if (purpose == "established" &
          focus == "structure" & n <= 200) {
        #C3: SRMR flex
        index <- "srmr"
        ap <- c(.001, .01, .05, .10)
        #fc <-
        #  suppressWarnings(sapply(ap, flex_co, fits = fits, index = index))
        fc <- suppressWarnings(lapply(
          ap, flex_co, fits = fits, index = index
        ))
        co <- vector("list", length(fc))
        gof <- co
        for (i in seq_len(length(fc))) {
          co[[i]] <- fc[[i]]$cutoff
          gof[[i]] <- fc[[i]]$gof
        }
        na <- fc[[1]]$`number of non-converging models`
        sh.na <- fc[[1]]$`share of non-converging models`
        #co <- fc[1,]
        #gof <- fc[4,]
        #na <- fc[6, 1]
        #sh.na <- fc[7, 1]
        rf <- fm[index]
      }
      if (purpose == "established" &
          focus == "structure" & n > 200) {
        #C4: SRMR flex
        index <- "srmr"
        ap <- c(.001, .01, .05, .10)
        #fc <-
        #  suppressWarnings(sapply(ap, flex_co, fits = fits, index = index))
        fc <- suppressWarnings(lapply(
          ap, flex_co, fits = fits, index = index
        ))
        co <- vector("list", length(fc))
        gof <- co
        for (i in seq_len(length(fc))) {
          co[[i]] <- fc[[i]]$cutoff
          gof[[i]] <- fc[[i]]$gof
        }
        na <- fc[[1]]$`number of non-converging models`
        sh.na <- fc[[1]]$`share of non-converging models`
        #co <- fc[1,]
        #gof <- fc[4,]
        #na <- fc[6, 1]
        #sh.na <- fc[7, 1]
        rf <- fm[index]
      }
      if (purpose == "novel" & focus == "cfa" & n <= 200) {
        #C5: CFI fix
        index <- "cfi"
        ap <- NULL
        fc <- NULL
        na <- NULL
        sh.na <- NULL
        co <- .95
        gof <- index_guess(index)
        rf <- fm[index]
      }
      if (purpose == "novel" & focus == "cfa" & n > 200) {
        #C6: SRMR flex
        index <- "srmr"
        ap <- c(.001, .01, .05, .10)
        #fc <-
        #  suppressWarnings(sapply(ap, flex_co, fits = fits, index = index))
        fc <- suppressWarnings(lapply(
          ap, flex_co, fits = fits, index = index
        ))
        co <- vector("list", length(fc))
        gof <- co
        for (i in seq_len(length(fc))) {
          co[[i]] <- fc[[i]]$cutoff
          gof[[i]] <- fc[[i]]$gof
        }
        na <- fc[[1]]$`number of non-converging models`
        sh.na <- fc[[1]]$`share of non-converging models`
        #co <- fc[1,]
        #gof <- fc[4,]
        #na <- fc[6, 1]
        #sh.na <- fc[7, 1]
        rf <- fm[index]
      }
      if (purpose == "novel" & focus == "structure" & n <= 200) {
        #C7: CFI & SRMR fix
        index <- c("cfi", "srmr")
        ap <- NULL
        fc <- NULL
        na <- NULL
        sh.na <- NULL
        co <- c(.95, .09)
        #gof <- sapply(index, index_guess)
        gof <- vapply(index, index_guess, character(1))
        rf <- fm[index]
      }
      if (purpose == "novel" & focus == "structure" & n > 200) {
        #C8: SRMR flex
        index <- "srmr"
        ap <- c(.001, .01, .05, .10)
        #fc <-
        #  suppressWarnings(sapply(ap, flex_co, fits = fits, index = index))
        fc <- suppressWarnings(lapply(
          ap, flex_co, fits = fits, index = index
        ))
        co <- vector("list", length(fc))
        gof <- co
        for (i in seq_len(length(fc))) {
          co[[i]] <- fc[[i]]$cutoff
          gof[[i]] <- fc[[i]]$gof
        }
        na <- fc[[1]]$`number of non-converging models`
        sh.na <- fc[[1]]$`share of non-converging models`
        #co <- fc[1,]
        #gof <- fc[4,]
        #na <- fc[6, 1]
        #sh.na <- fc[7, 1]
        rf <- fm[index]
      }
    }
    if (override) {
      if (is.null(index))
        stop(
          "You chose to override the recommendation patterns. In this case, you need to provide index names. Please revise."
        )
      index <- tolower(index)
      ap <- c(.001, .01, .05, .10)
      #fc <-
      #  suppressWarnings(sapply(ap, flex_co, fits = fits, index = index))
      fc <- suppressWarnings(lapply(
        ap, flex_co, fits = fits, index = index
      ))
      co <- vector("list", length(fc))
      gof <- co
      for (i in seq_len(length(fc))) {
        co[[i]] <- fc[[i]]$cutoff
        gof[[i]] <- fc[[i]]$gof
      }
      na <- fc[[1]]$`number of non-converging models`
      sh.na <- fc[[1]]$`share of non-converging models`
      #co <- fc[1,]
      #gof <- fc[4,]
      #na <- fc[6, 1]
      #sh.na <- fc[7, 1]
      rf <- fm[index]
    }
    if (!is.null(ap)) {
      tab <-
        data.frame(matrix(
          unlist(co),
          nrow = length(ap),
          ncol = length(index),
          byrow = TRUE
        ))
      names(tab) <- toupper(index)
      rownames(tab) <- ap
      gof <-
        data.frame(matrix(
          unlist(gof),
          nrow = length(ap),
          ncol = length(index),
          byrow = TRUE
        ))
      names(gof) <- toupper(index)
      rownames(gof) <- ap
      decs <- gof
      for (i in seq_len(length(index))) {
        decs[, i] <- rf[i] - tab[, i]
        decs[, i] <-
          ifelse(gof[, i] == TRUE, decs[, i], decs[, i] * -1)
        decs[, i] <- ifelse(decs[, i] > 0, "confirmed", "rejected")
      }
      rownames(decs) <- paste0("cutoff ", ap)
      rownames(tab) <- paste0("cutoff ", ap)
      #gof <- data.frame(unname(sapply(index, index_guess)))
      gof <- data.frame(unname(vapply(index, index_guess, character(1))))
      gof <- cbind(gof, unname(rf))
      rownames(gof) <- toupper(index)
      names(gof) <- c("type", "fit.values")
      gof[, "fit.values"] <-
        round(gof[, "fit.values"], digits = digits)
      tab <- round(tab, digits = digits)
      if (fits$rep < 500)
        warning(
          "The number of replications is lower than the recommended minimum of 500. Consider with care."
        )
      res <- list(
        "recommended" = gof,
        "cutoffs" = tab,
        "decisions" = decs,
        "replications" = fits$rep,
        "comment" = ifelse(
          override,
          "Override mode",
          "Recommendations based on flexible cutoffs and Mai et al. (2021)"
        )
      )
    }
    if (is.null(ap)) {
      decs <- rf - co
      decs <- ifelse(gof == "GoF", decs, decs * -1)
      decs <- ifelse(decs > 0, "confirmed", "rejected")
      #gof <- data.frame(unname(sapply(index, index_guess)))
      gof <- data.frame(unname(vapply(index, index_guess, character(1))))
      gof <- cbind(gof, unname(rf))
      rownames(gof) <- toupper(index)
      names(gof) <- c("type", "fit.values")
      gof[, "fit.values"] <-
        round(gof[, "fit.values"], digits = digits)
      res <- list(
        "recommended" = gof,
        "decisions" = decs,
        "comment" = "Recommendations based on fixed cutoffs and Mai et al. (2021)"
      )
    }
    return(res)
  }
