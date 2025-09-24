#' Obtain fit statistics from correctly and misspecified models.

#' @details
#' Originally proposed by Niemand & Mai (2018), flexible cutoffs (hereafter FCO1) have been developed with only a correct model in mind. This simple, first approach of simulated cutoffs for fit indices hence cannot consider misspecified models under consideration. To improve on this, multiple decision rules are introduced that – based on a misspecified model – allow to determine the Type II error of misfit for a given cutoff and fit indicator. With this feature in mind, this enables different decision rules (i.e., what type of error is considered and how) based on pertinent literature:
#' FCO1: The original proposed flexible cutoffs (Niemand & Mai, 2018), only considering Type I error. The p (e.g., 5%) quantile of all simulated correct models is taken (assuming a GoF like CFI, alternatively 1-p for a BoF like SRMR).
#' FCO2: A modified flexible cutoff considering Type I and II error. The 1-p (e.g., 95%) quantile of all simulated misfit models is taken when this quantile is smaller than or equal to the p (e.g., 5%) quantile. Otherwise, the p (e.g., 5%) quantile of all simulated correct models is taken. FCO2 (alike FCO1) always provides a cutoff (assuming a GoF like CFI, alternatively p (misfit model) and 1-p (correct model) for a BoF like SRMR)
#' DFI: A modified dynamic cutoff considering Type I and II error (McNeish & Wolff, 2023). The 1-p (e.g., 95%) quantile of all simulated misfit models is taken when this quantile is smaller than or equal to the p (e.g., 5%) quantile. Otherwise, no cutoff is provided (NA). The DFI-decision rule tends to provide no cutoff when correct and misfit models overlap strongly (assuming a GoF like CFI, alternatively p (misfit model) and 1-p (correct model) for a BoF like SRMR).
#' CP: A modified cutoff considering Type I and II error via an optimal cutpoint (Groskurth et al., 2022). The cutoff is found by taking the cutoff with the highest sum of sensitivity (1 – Type I error) and 1 – specificity (Type II error) derived from simulated correct and misfit models via cutpointr::cutpointr. CP always provides a cutoff.
#' Fix: Fixed cutoffs (Hu & Bentler, 1999) are also provided for comparison.
#' Since it cannot be objectively determined what level or type of misspecification (and to which extent) demarcates "acceptable" and "unacceptable" misfit, generalizability of the misspecification procedure becomes a vital question. To overcome this issue, a PROCESS-like (Hayes, 2017) approach is proposed. Instead of expecting that the user provides an appropriate misspecified model (as in simsem), which might be highly user-unfriendly, the user only provides a type of misfit model via model.type. This argument defines the number of structural model misspecifications (first integer), measurement model misspecifications (second integer) and residual covariance misspecifications (third integer) assumed for the misfit model. For example, "100" refers to one structural model misspecification (factor correlations set to 0), zero measurement model misspecifications (no cross-loadings set to 0) and zero residual covariances introduced to the correct model described in mod. The default is set to "111", which corresponds to a model where one correlation, one cross-loading and one residual covariance may be overlooked. If a researcher is certain, that only one type of misspecification is important, the value can be changed to a "100", "010", or "001" model for example. Comparing multiple model.type specifications is recommended.
#' Based on feedback on the previous versions of FCO, some new features are integrated, most importantly direct input of a fitted lavaan object (fit), a one-step calculation of fits, support for different types of variables, extensive checking of the model and data characteristics by an internal function providing better warning and error descriptions, directly setting skewness (sk) and kurtosis (ku), easier parallelization options and random generation (random) of the misfit model. Further, two more functions have been introduced to better compare the implications (e.g., in terms of implied Type I and II errors) across decision rules and to visualize the simulated correct and misfit model distributions.
#' Function flex_co2 now provides detailed information on the quantiles, cutoffs (all decision rules plus fixed cutoffs), evaluation (including sum of errors, Type I error, Type II error estimates from the simulation data), a notation and overlap statistics (percentage overlap from overlapping::overlap, AUC from cutpointr::cutpointr and a U-test converted to d from rcompanion::wilcoxonR) to identify the degree of overlap between simulated correct and misfit model distributions. Thereby, users can select an appropriate cutoff, inspect which decision rule works best for their data and model and determine how much the distributions overlap. Function plot_fit2 complements these features and plots the simulated distributions for correct and misfit models to illustrate the distributions, their overlap and the cutoff quantiles.
#' Multiple notes are given regarding the estimates 1) "DFI" does not refer to using the same approach as dynamic::cfaHB as proposed by McNeish & Wolf (2023). DFI only refers to the authors' decision rule, but applying the same PROCESS-like approach as FCO1, FCO2 and CP. 2). For compatibility reasons, functions gen_fit, flex_co, pop_mod, recommend and recommend_dv are (virtually) unchanged and only apply for FCO1. Likewise, functions flex_co2 and plot_fit2 only apply for results of this function gen_fit2.
#' @param fit A fitted object from lavaan.
#' @param mod A lavaan model to specify the CFA. If a fitted lavaan object is provided, the model is taken from this object via lavaan::parTable.
#' @param x A dataset for the model of nrow observations and ncol indicators. If a fitted lavaan object is provided, the dataset is taken from this object.
#' @param n A sample size specified instead of a dataset. Requires a population model via pop.mod. If a fitted lavaan object is provided, the sample size is taken from this object.
#' @param model.type A model type defining the number of structural model misspecifications (first integer), measurement model misspecifications (second integer) and residual covariance misspecifications (third integer) assumed for the misspecified model. Can be written as a character (default: "111") or as an integer vector (default: c(1,1,1)). For each integer, a maximum of 3 is supported so far. For example, "100" refers to one structural model misspecification (factor correlations set to 0), zero measurement model misspecifications (cross-loadings set to 0) and zero residual covariances introduced to the correct model described in mod. A "000" model is not allowed as it would be identical to the correct model.
#' @param pop.mod For flexibility reasons, an optional lavaan population model can be provided. This is required together with n if x is missing or no fitted lavaan object is used.
#' @param rep Number of replications to be simulated (default: 500).
#' @param type Type of underlying population model. Based on the model provided, a population model is derived to simulate the fit indices by function the internal function pop_mod. The type determines the factor loadings and covariances assumed for this population model. NM (the default): Uses the factor loadings and covariances from Niemand & Mai's (2018) simulation study. HB: Uses the factor loadings and covariances from Hu & Bentler's (1999) simulation study. EM: Empirical uses the given factor loadings and covariances. The underlying population model is used to derive a misspecified population model based on the model provided in model.type.
#' @param cfa If TRUE (the default), the population model is generated assuming a CFA, i.e., there are only loadings and correlations. If FALSE, the model can be a regression model or any other type of SEM. In this case, the argument type determines how the population model is built. In case of “EM”, the population model is defined from the model with the given parameters. If type is “NM” or “HB”, the population model is defined based on the loadings, correlations and regression coefficients (beta) given effect size values for these three (loadings: see es.lam; correlations: see es.cor; beta: see es.f2). In the latter cases, aco and afl are not used.
#' @param data.types Types of the manifest variables. Users can specify a vector of the length of variables with C = count, B = binary, O = ordinal, N = normal in the same order as in the dataset. These types are then used to simulate data for the population model based on median (count), mean (binary), cumulative relative frequencies of values (ordinal) as well as mean and variance (normal) applying the PoisBinOrdNor::intermat and PoisBinOrdNor::genPBONdata functions from package PoisBinOrdNor. That is, categorial and binary variables are also supported. Argument type is set to EM when data.types are defined (otherwise, normal data would be implied).
#' @param esti The estimator to be used for model estimation in lavaan, defaults to "ML". Consider changing when needed. If a fitted lavaan object is provided, the estimator is taken from this object.
#' @param cores How many cores should be used for multiple cores? The default is 2. Consider the available number of cores of your system.
#' @param standardized Are factor loadings assumed to be standardized and covariances to be correlations (default: NULL)? The internal function pop_mod checks this feature and returns a warning if set to TRUE (any > 1) or FALSE (all < 1). Otherwise, TRUE or FALSE is guessed from the loadings (all < 1 leads to TRUE, any > 1 to FALSE).
#' @param es.lam Effect size assumed for the loadings when cfa = FALSE and type is not “EM”. Options are ‘low’ (.7), ‘moderate’ (.8) and ‘large’ (.9). Defaults to ‘low’.  The loadings are equal for type “NM” and vary according to Hu & Bentler (1999) for “HB”.
#' @param es.cor Effect size assumed for the correlations when cfa = FALSE and type is not “EM”. Options are ‘low’ (.1), ‘moderate’ (.3) and ‘large’ (.5) based on Cohen’s (1988) conventions. Defaults to ‘large’.  The correlations are equal for type “NM” and vary according to Hu & Bentler (1999) for “HB”.
#' @param es.f2 Effect size assumed for the regression coefficients when cfa = FALSE and type is not “EM”. Options are ‘low’ (.02), ‘moderate’ (.15) and ‘large’ (.35) based on Cohen’s (1988) conventions. Defaults to ‘large’.  The regression coefficients are equal for type “NM” and “HB”.
#' @param sk Should skewness (default: 0) be assumed to indicate multivariate normality? In case of excessive skewness via psych::skew, a warning is provided and the user can enter a more appropriate sk for the data in the next run.
#' @param ku Should kurtosis (default: 1) be assumed to indicate multivariate normality? In case of excessive kurtosis via psych::kurtosi, a warning is provided and the user can enter a more appropriate ku for the data in the next run.
#' @param seed The seed to be set to obtain reproducible cutoffs (default: 1111). Defines a vector of length rep with the seed being the first value.
#' @param random Should the misspecified population model be generated randomly (default: TRUE)? To avoid a bias by always misspecifying the same parameter in the model based on model.type for every replication (= FALSE), the parameter can be randomly selected (= TRUE). Results differ slightly, yet random = TRUE is a bit slower.
#' @return A list of simulated fit statistics for correct models (correct.fits) and misfit models (miss.fits).
#' @references Groskurth, K., Bhaktha, N., & Lechner, C. (2022). Making model judgments ROC (K)-solid: Tailored cutoffs for fit indices through simulation and ROC analysis in structural equation modeling. https://psyarxiv.com/62j89/download?format=pdf
#' @references Hu, L., & Bentler, P. M. (1999). Cutoff criteria for fit indexes in covariance structure analysis: Conventional criteria versus new alternatives. Structural Equation Modeling, 6(1), 1–55. https://doi.org/10.1080/10705519909540118
#' @references McNeish, D., & Wolf, M. G. (2023). Dynamic fit index cutoffs for confirmatory factor analysis models. Psychological Methods, 28(1), 61–88. https://doi.org/10.1037/met0000425
#' @references Niemand, T., & Mai, R. (2018). Flexible cutoff values for fit indices in the evaluation of structural equation models. Journal of the Academy of Marketing Science, 46(6), 1148–1172. https://doi.org/10.1007/s11747-018-0602-9
#' @importFrom stats cor cov median na.omit quantile var
#' @importFrom dplyr %>% as_tibble
#' @importFrom parallel mcmapply
#' @importFrom utils combn
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
#' flex_co2(fits)
#' plot_fit2(fits)
#'
#' #Different data types
#' dat <- lavaan::HolzingerSwineford1939
#' cdat <- dplyr::select(dat, x1, x2, x3, x4, x5, x6, x7, x8, x9)
#' #For demo purposes, some variables are changed:
#' cdat <- cdat %>% mutate(
#'   x1 = round(x1, digits = 0),
#'   x3 = round(x3, digits = 0),
#'   x2 = ifelse(x2 > 4, 1, 0),
#'   x4 = ifelse(x4 > 2, 1, 0),
#'   x5 = round(x5, digits = 0),
#'   x8 = round(x8, digits = 0)
#' )
#' cfit <- cfa(model = HS.model, data = cdat)
#' #Note: Demonstration only! Please use higher numbers of replications for your applications (>= 500).
#' cfits <- gen_fit2(fit = cfit,
#'                   data.types = c("C", "B", "C", "B", "O", "N", "N", "O", "N"), rep = 10)
#' flex_co2(cfits)
#' plot_fit2(cfits)
#'
#' #Multiple fit indices
#' flex_co2(fits, index = c("cfi", "SRMR", "RMSEA"))
#' plot_fit2(fits, index = c("cfi", "SRMR", "RMSEA"))
#' @export
gen_fit2 <- function(fit = NULL,
                     mod = NULL,
                     x = NULL,
                     model.type = "111",
                     pop.mod = NULL,
                     n = NULL,
                     rep = 500L,
                     type = "NM",
                     cfa = TRUE,
                     data.types = NULL,
                     esti = "ML",
                     cores = 2,
                     standardized = TRUE,
                     es.lam = "low",
                     es.cor = "large",
                     es.f2 = "moderate",
                     sk = 0,
                     ku = 1,
                     seed = 1111,
                     random = TRUE) {
  #Get results from fit object
  if (!is.null(fit)) {
    lvo <- get_lvo(fit)
    mod <- lvo$mod
    n <- lvo$n
    esti <- lvo$esti
    x <- lvo$x
  }

  #Compatibility with NM
  if (type == "NM") {
    afl <- .7
    aco <- .3
  }

  #Model handling
  if (is.null(fit)) {
    if (is.null(pop.mod) && is.null(mod))
      stop("Provide model (mod) or pop.mod or fit.")
    if (is.null(x) && is.null(n))
      stop("Provide data (x), sample size (n) or fit.")
    if (is.null(n))
      n <- nrow(x)
  }

  #Get par.type - omitted
  par.type <- "foreach"

  #Seeding
  if (is.null(seed))
    seed <- sample(1111:9999, 1L)
  RNGkind("L'Ecuyer-CMRG")
  seeds <- seed:(seed + rep - 1L)

  #Get pop.mod
  if (!is.null(pop.mod)) {
    cpop <- pop.mod
    if (is.null(mod))
      mod <- get_free_mod(cpop)
  }

  if (is.null(pop.mod)) {
    if (cfa) {
      cpop <- pop_mod(
        mod = mod,
        x = x,
        type = type,
        standardized = standardized,
        afl = afl,
        aco = aco,
        data.types = data.types,
        seed = seed
      )$pop.mod
    } else {
      if (!is.null(data.types))
        x <- gen_nnd(x, data.types = data.types, seed = seed)
      cpop <- pop_mod_reg(
        mod = mod,
        x = x,
        type = type,
        es.lam = es.lam,
        es.cor = es.cor,
        es.f2 = es.f2
      )$pop.mod
    }
  }

  #Generate data if needed
  if (is.null(x)) {
    x <- lavaan::simulateData(
      model = pop.mod,
      model.type = "cfa",
      std.lv = TRUE,
      auto.fix.first = FALSE,
      sample.nobs = n,
      seed = seed,
      kurtosis = ku,
      skewness = sk
    )
  }

  #Miss model building
  if (is.character(model.type))
  model.type <- unname(as.numeric(strsplit(model.type, "")[[1]]))
  sm <- model.type[1]
  mm <- model.type[2]
  rcv <- model.type[3]

  #Check
  checker(
    pop.mod = cpop,
    cfa = cfa,
    x = x,
    esti = esti,
    sm = sm,
    mm = mm,
    rcv = rcv,
    cores = cores
  )

  #Sim correct
  cfit <- run_parallel(
    mini.generator,
    seeds,
    cores,
    par.type,
    pop.mod = cpop,
    mod = mod,
    n = n,
    sk = sk,
    ku = ku,
    esti = esti
  )

  #Sim miss
  if (!random) {
    mpop <- build_and_check(seed, sm, mm, rcv, cpop, cfa, x, esti, cores)
    mfit <- run_parallel(
      mini.generator,
      seeds,
      cores,
      par.type,
      pop.mod = mpop,
      mod = mod,
      n = n,
      sk = sk,
      ku = ku,
      esti = esti
    )
  } else {
    gen_miss <- function(seed) {
      mpop <- build_and_check(seed, sm, mm, rcv, cpop, cfa, x, esti, cores)
      mini.generator(
        pop.mod = mpop,
        seed = seed,
        mod = mod,
        n = n,
        sk = sk,
        ku = ku,
        esti = esti
      )
    }
    mfit <- run_parallel(gen_miss, seeds, cores, par.type)
  }

  #Output
  return(list(
    correct.fits = dplyr::as_tibble(cfit, .name_repair = "minimal"),
    miss.fits = dplyr::as_tibble(mfit, .name_repair = "minimal")
  ))
}
