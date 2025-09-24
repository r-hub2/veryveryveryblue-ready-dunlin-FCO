#' Obtain fit statistics from one or two models
#'
#' @param mod1 A lavaan model to specify the CFA.
#' @param mod2 Another lavaan model for a model comparison. If missing and merge.mod = TRUE, a merged model from function merge_factors is estimated based on mod1.
#' @param x A dataset for the model of nrow observations (minimum: 50) and ncol indicators (minimum: 4)
#' @param n A sample size specified instead of a dataset (minimum: 50, maximum: 50000). Requires a population model via pop.mod1.
#' @param rep Number of replications to be simulated (default: 500, minimum: 10, maximum: 5000)
#' @param type Type of underlying population model. Based on the model(s) provided, a population model is derived to simulate the fit indices by function pop_mod. The type determines
#' the factor loadings and covariances assumed for this population model. NM (the default when only one model is provided): Uses the factor loadings and
#' covariances from Niemand & Mai's (2018) simulation study. HB: Uses the factor loadings and covariances from Hu & Bentler's (1999) simulation study.
#' EM: Empirical (the default when two models are provided or merge.mod is TRUE), uses the given factor loadings and covariances.
#' @param dv Should the fit statistics be calculated for discriminant validity testing? If no (the default), this is not assumed. If yes, consider the arguments of merge.mod, dv.factors and cutoff.
#' So far, two options of discriminant validity testing are supported. Constraining: A factor correlation between two factors can be constrained as selected by the dv.factors argument. In this case, dv.cutoff applies and merge.mod is not required.
#' Merging: Two factors can be merged into one, again controlled by the dv.factors argument. In this case, merge.mod applies and dv.cutoff is not required (as cutoff = 1 is implied).
#' @param dv.factors Names of the factors to be considered. Must be equal to 2. If missing (the default), the first and second factor of the model are selected.
#' @param merge.mod This is used for merging. If FALSE (the default), fit measures for mod1 are estimated for a single model as long as no mod2 is provided. If TRUE, a merged model from function merge_factors is estimated based on mod1. In this case, no mod2 is required.
#' @param dv.cutoff This is used for constraining. It determines the critical correlation assumed to be a cutoff for discriminant validity testing.
#' For example, based on Rönkkö & Cho (2020), a cutoff of .9 indicates a severe issue in discriminant validity between the selected factors. Cutoffs between .8 and 1 are recommended.
#' The function returns a warning, if the cutoff is below .8.
#' @param standardized Are factor loadings assumed to be standardized and covariances to be correlations (default: TRUE)?
#' @param assume.mvn Should multivariate normality (mvn) be assumed? If TRUE (the default), kurtosis and skewness are set to 1 for simulated data.
#' If FALSE, kurtosis and skewness are estimated from dataset x via semTools::mardiaKurtosis and semTools::mardiaSkew.
#' @param multi.core Should multiple cores be used to simulate fit indices?
#' If TRUE (the default), mclapply (on Linux or Mac machines) or parLapply (on Windows machines) from parallel package with the number of specified cores is used. If FALSE, a single core is used.
#' @param cores How many cores should be used for multiple cores? The default is 2. Consider the available number of cores of your system.
#' @param seed The seed to be set to obtain reproducible cutoffs (default: 1111). Defines a vector of length rep with the seed being the first value.
#' @param pop.mod1 For flexibility reasons, an optional lavaan population model can be provided. This is required together with n if x is missing.
#' @param pop.mod2 Another optional lavaan population model.
#' @return A list of simulated fit statistics (fco) and all previously defined parameters.
#' @references Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure analysis: Conventional criteria versus new alternatives. Structural Equation Modeling, 6(1), 1–55. https://doi.org/10.1080/10705519909540118
#' @references Niemand, T., & Mai, R. (2018). Flexible cutoff values for fit indices in the evaluation of structural equation models. Journal of the Academy of Marketing Science, 46(6), 1148—1172. https://doi.org/10.1007/s11747-018-0602-9
#' @references Rönkkö, M., & Cho, E. (2020). An updated guideline for assessing discriminant validity. Organizational Research Methods. https://doi.org/10.1177/1094428120968614
#' @examples
#'#Note: Demonstration only! Please use higher numbers of replications for your applications (>= 500).
#' #A single model to obtain fit indices for
#'mod <- "
#'F1 =~ Q5 + Q7 + Q8
#'F2 =~ Q2 + Q4
#'F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
#'F4 =~ Q1 + Q17
#'F5 =~ Q6 + Q14 + Q15 + Q16
#'"
#' fits.single <- gen_fit(mod1 = mod, x = bb1992, rep = 10, standardized = FALSE)
#'
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

#' #Two models for discriminant validity testing, this resembles constraining with a cutoff of .9
#'fits.dv.con <- gen_fit(
#'  mod1 = mod,
#'  x = bb1992,
#'  rep = 10,
#'  dv = TRUE,
#'  dv.factors = c("F4", "F5"),
#'  dv.cutoff = .9
#')
#'
#' #Two models for discriminant validity testing, this resembles merging.
#'fits.dv.merge <- gen_fit(
#'  mod1 = mod,
#'  x = bb1992,
#'  rep = 10,
#'  dv = TRUE,
#'  dv.factors = c("F4", "F5"),
#'  merge.mod = TRUE
#')
#'}
#' @export
gen_fit <-
  function(mod1 = NULL,
           mod2 = NULL,
           x = NULL,
           n = NULL,
           rep = 500,
           type = "NM",
           dv = FALSE,
           dv.factors = NULL,
           merge.mod = FALSE,
           dv.cutoff = .9,
           standardized = TRUE,
           assume.mvn = TRUE,
           multi.core = TRUE,
           cores = 2,
           seed = 1111,
           pop.mod1 = NULL,
           pop.mod2 = NULL) {
    #Checks
    checkmate::assertCharacter(mod1,
                               fixed = "=~", null.ok = TRUE)
    checkmate::assertCharacter(mod2,
                               fixed = "=~", null.ok = TRUE)
    checkmate::assertDataFrame(x,
                               min.rows = 50,
                               min.cols = 4,
                               col.names = "unique",
                               null.ok = TRUE)
    checkmate::assertNumeric(n, lower = 50, upper = 50000, null.ok = TRUE)
    checkmate::assertNumeric(rep, lower = 10, upper = 5000)
    checkmate::assert(
      checkmate::checkCharacter(type,
                                pattern = "EM"),
      checkmate::checkCharacter(type,
                                pattern = "HB"),
      checkmate::checkCharacter(type,
                                pattern = "NM")
    )
    checkmate::assertLogical(dv)
    checkmate::assertVector(dv.factors, len = 2, null.ok = TRUE)
    checkmate::assertLogical(merge.mod)
    checkmate::assertVector(dv.cutoff, len = 1, null.ok = FALSE)
    if (dv.cutoff > 1 |
        dv.cutoff < 0)
      stop("Cutoff not between 0 and 1. Please revise.")
    if (dv.cutoff < .8)
      warning("Cutoff below recommended standards for discriminant validity. Consider with care.")
    checkmate::assertLogical(standardized)
    checkmate::assertLogical(assume.mvn)
    checkmate::assertLogical(multi.core)
    checkmate::assertNumeric(cores, lower = 1)
    checkmate::assertNumeric(seed, any.missing = FALSE)
    #Check the structure of pop.mods
    if (!is.null(pop.mod1) &
        is.list(pop.mod1))
      pop.mod1 <- pop.mod1$pop.mod
    if (!is.null(pop.mod2) &
        is.list(pop.mod2))
      pop.mod2 <- pop.mod2$pop.mod
    checkmate::assertCharacter(pop.mod1,
                               fixed = "=~",
                               null.ok = TRUE)
    checkmate::assertCharacter(pop.mod2,
                               fixed = "=~",
                               null.ok = TRUE)
    stopn <- "Please either provide a) model and dataset or b) population model and sample size."
    if (is.null(x)) {
      if (is.null(n)) stop(stopn)
      if (!is.null(n)) {
        if(is.null(pop.mod1)) stop(stopn)
      }
    }
    if (!is.null(x)) {
      if(is.null(mod1)) stop(stopn)
      if(!is.null(pop.mod1)) stop(stopn)
    }
    if (!is.null(pop.mod1)) {
      if (!grepl("*", pop.mod1, fixed = TRUE)) stop("You provided a population model that looks like a model. Please revise.")
    }
    if (!is.null(pop.mod2)) {
      if (!grepl("*", pop.mod2, fixed = TRUE)) stop("You provided a second population model that looks like a model. Please revise.")
    }
        #Check the current no. of fit indices provided by lavaan
    # nf <-
    #   length(lavaan::fitmeasures(
    #     lavaan::cfa(
    #       mod = "F =~ x1 +x2 +x3",
    #       data = lavaan::HolzingerSwineford1939,
    #       estimator = "MLM"
    #     )
    #   ))
    nf <- 69
    #Set normality
    if (!assume.mvn & !is.null(x)) {
      k <-
        unname(semTools::mardiaKurtosis(x)["b2d"]) / (ncol(x) * (ncol(x) + 2))
      #Centered kurtosis
      s <- unname(semTools::mardiaSkew(x)["b1d"])
    }
    if (assume.mvn) {
      k <- 1
      s <- 1
    }
    seeds <- seq(from = seed, to = seed + rep - 1)
    RNGkind("L'Ecuyer-CMRG")
    #Seeds are hold constant for replication
    if (multi.core) {
      if (parallel::detectCores() < cores)
        stop(
          "The number of available cores is lower than the number of selected cores. Please revise."
        )
      if (.Platform$OS.type == "windows") {
        clus <- parallel::makePSOCKcluster(cores)
      }
    }
    if (!multi.core) {
      cores <- 1
    }
    g <- expand.grid(nullm2 = c(TRUE, FALSE),
                     mm = c(TRUE, FALSE),
                     dv = c(TRUE, FALSE))
    g$opt <- c(rep("dv", 5), rep("cfa", 3))
    g$mode <-
      c(
        "merging",
        "constraining",
        "constraining",
        "constraining",
        "merging",
        "dual",
        "single",
        "dual"
      )
    wh <- which(g$nullm2 == is.null(mod2) &
                  g$mm == merge.mod & g$dv == dv)
    if (wh == 2 | wh == 3 | wh == 4)
      message(
        "Constraining is selected as the discriminant validity testing option given the provided arguments."
      )
    if (wh == 1)
      message(
        "Merging is selected as the discriminant validity testing option given the provided arguments."
      )
    if (wh == 5)
      message("Merge.factors is TRUE, but dv is FALSE. Merging is selected while assuming dv to be TRUE.")
    if (wh == 6)
      message(
        "Two models are provided, but merge.factors is TRUE. Estimation continues while assuming merge.factors to be FALSE."
      )
    mode <- g[wh, "mode"]
    if (mode == "single") {
      if (is.null(pop.mod1)) {
        pop.mod <-
          pop_mod(
            mod = mod1,
            x = x,
            type = type,
            standardized = standardized
          )$pop.mod
      }
      if (!is.null(pop.mod1)) pop.mod <- pop.mod1
      if (identical(mod1, pop.mod)) stop("You provided a model that looks like a population model. Please revise.")
      free <- mod1
      if (multi.core & .Platform$OS.type != "windows") {
        fits <-
          parallel::mclapply(1:rep, function(r) {
            generator(
              x = x,
              n = n,
              seed = seeds[r],
              mode = mode,
              pop.mod1 = pop.mod,
              free1 = free,
              s = s,
              k = k,
              nf = nf
            )
          }, mc.cores = cores)
      }
      if (multi.core & .Platform$OS.type == "windows") {
        fits <-
          parallel::parLapply(clus, 1:rep, function(r) {
            generator(
              x = x,
              n = n,
              seed = seeds[r],
              mode = mode,
              pop.mod1 = pop.mod,
              free1 = free,
              s = s,
              k = k,
              nf = nf
            )
          })
      }
      if (!multi.core) {
        fits <-
          lapply(1:rep, function(r) {
            generator(
              x = x,
              n = n,
              seed = seeds[r],
              mode = mode,
              pop.mod1 = pop.mod,
              free1 = free,
              s = s,
              k = k,
              nf = nf
            )
          })
      }
    }
    if (mode != "single") {
      type <- "EM"
      if (is.null(pop.mod1)) {
        pop.mod1 <-
          pop_mod(
            mod = mod1,
            x = x,
            type = type,
            standardized = standardized
          )$pop.mod
      }
      if (identical(mod1, pop.mod1)) stop("You provided a model that looks like a population model. Please revise.")
      free1 <- mod1
      if (mode == "dual" | mode == "constraining") {
        #Dual models, no dv
        if (is.null(pop.mod2) & mode == "dual") {
          if (is.null(mod2))
            stop("No second model provided. Please revise.")
          free2 <- mod2
        }
        #Constraining:
        if (is.null(pop.mod2) & mode == "constraining") {
          free2 <-
            ifelse(
              is.null(mod2),
              constr_mod(mod1, dv.factors = dv.factors, dv.cutoff = dv.cutoff),
              mod2
            )
          mod2 <- free2
        }
        if (!is.null(pop.mod2) & mode == "constraining") {
          #This is required as otherwise the simulation would take fit measures from the population model without the cutoff.
          free2 <-
            ifelse(is.null(mod2),
                   get_free(lavaan::lavaanify(
                     pop_mod_dv(
                       pop.mod2,
                       dv.factors = dv.factors,
                       dv.cutoff = dv.cutoff
                     )
                   ), dv.factors, mode),
                   mod2)
          mod2 <- free2
        }
      }
      #Merging:
      if (mode == "merging") {
        pt2 <-
          try(merge_factors(lavaan::cfa(pop.mod1, x, warn = FALSE)), silent = TRUE)
        if (inherits(pt2, "try-error"))
          stop("Merging two-factors not successful. Please check.")
        free2 <-
          ifelse(is.null(mod2), get_free(pt2, dv.factors, mode), mod2)
        mod2 <- free2
      }
      if (multi.core & .Platform$OS.type != "windows") {
        fits <-
          parallel::mclapply(1:rep, function(r) {
            generator(
              x = x,
              n = n,
              seed = seeds[r],
              mode = mode,
              pop.mod1 = pop.mod1,
              free1 = free1,
              free2 = free2,
              s = s,
              k = k,
              nf = nf
            )
          }, mc.cores = cores)
      }
      if (multi.core & .Platform$OS.type == "windows") {
        fits <-
          parallel::parLapply(clus, 1:rep, function(r) {
            generator(
              x = x,
              n = n,
              seed = seeds[r],
              mode = mode,
              pop.mod1 = pop.mod1,
              free1 = free1,
              free2 = free2,
              s = s,
              k = k,
              nf = nf
            )
          })
      }
      if (!multi.core) {
        fits <-
          lapply(1:rep, function(r) {
            generator(
              x = x,
              n = n,
              seed = seeds[r],
              mode = mode,
              pop.mod1 = pop.mod1,
              free1 = free1,
              free2 = free2,
              s = s,
              k = k,
              nf = nf
            )
          })
      }
    }
    if (multi.core & .Platform$OS.type == "windows") {
      parallel::stopCluster(clus)
    }
    res <- list(
      "fco" = fits,
      "mod1" = mod1,
      "mod2" = mod2,
      "x" = x,
      "n" = n,
      "rep" = rep,
      "type" = type,
      "dv" = dv,
      "dv.factors" = dv.factors,
      "merge.mod" = merge.mod,
      "dv.cutoff" = dv.cutoff,
      "standardized" = standardized,
      "assume.mvn" = assume.mvn,
      "multi.core" = multi.core,
      "cores" = cores,
      "seed" = seed,
      "pop.mod1" = pop.mod1,
      "pop.mod2" = pop.mod2
    )
    return(res)
  }
