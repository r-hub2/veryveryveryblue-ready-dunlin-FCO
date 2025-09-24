#' Plotting the distributions of selected simulated fit indices
#' @details For details, please refer to gen_fit2. Please note that the results are only based on the simulation of the model specified and its misspecified variant under the conditions of the gen_fit2 function.
#' @param fits The object returned from gen_fit2 which is a list of correct fit values and misspecified fit values.
#' @param correct.fits For compatibility reasons, the correct fit values can be defined separately.
#' @param miss.fits For compatibility reasons, the misspecified fit values can be defined separately.
#' @param index A vector of length >= 1 with names of fit indices the user wishes to explore. Capitalization does not matter, either e.g., CFI or cfi are accepted. Default is CFI as it might be easier to understand.
#' @param alpha The acceptable Type I error representing the empirical quantile p, see details in gen_fit2. Multiple values can be provided as a vector.  Default is .05.
#' @param beta The acceptable Type II error representing the empirical quantile p, see details in gen_fit2. Multiple values can be provided as a vector. Default is c(.05, .10).
#' @return A ggplot2 object with the simulated cutoffs for correct and misspecified models iterated across the number of indices provided.
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
#' #Default plot:
#' plot_fit2(fits)
#' #Changed alpha and beta values:
#' plot_fit2(fits, alpha = .05, beta = .05)
#' plot_fit2(fits, alpha = .10, beta = .20)
#' #Different fit indices:
#' plot_fit2(fits, index = c("CFI", "SRMR", "RMSEA"))
#' @export
plot_fit2 <-
  function(fits = NULL,
           correct.fits = NULL,
           miss.fits = NULL,
           index = "CFI",
           alpha = .05,
           beta = .10) {
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
    ind <- index
    tt <- vector("list", length(ind))
    for (i in 1:length(ind)) {
      tt[[i]] <-
        cbind(dplyr::select(cf, tolower(ind[i])), dplyr::select(mf, tolower(ind[i])))
      names(tt[[i]]) <- c("correct", "misspecified")
      tt[[i]] <-
        tt[[i]] %>% tidyr::pivot_longer(
          cols = dplyr::everything(),
          names_to = "Model",
          values_to = "fit"
        )
      tt[[i]]$index <- rep(ind[[i]], nrow(tt[[i]]))
    }
    gd <- do.call("rbind", tt)
    vdc <- gd %>%
      dplyr::filter(Model == "correct") %>%
      dplyr::reframe(qc = quantile(
        fit,
        probs = c(alpha, 1 - alpha),
        type = 8,
        na.rm = T
      ),
      .by = index)
    vdc$alpha <- rep(c(alpha, 1 - alpha), length(ind))
    vdc$type <- sapply(vdc$index, index_guess, simplify = T)
    vdc <-
      dplyr::bind_rows(
        vdc %>% dplyr::filter(alpha == {
          {
            alpha
          }
        } & type == "GoF"),
        vdc %>% dplyr::filter(alpha != {
          {
            alpha
          }
        } & type == "BoF")
      )

    vdm <- gd %>%
      dplyr::filter(Model == "misspecified") %>%
      dplyr::reframe(qm = quantile(
        fit,
        probs = c(beta, 1 - beta),
        type = 8,
        na.rm = T
      ),
      .by = index)
    vdm$beta <- rep(c(beta, 1 - beta), length(ind))
    vdm$type <- sapply(vdm$index, index_guess, simplify = T)
    vdm <-
      dplyr::bind_rows(
        vdm %>% dplyr::filter(beta == {
          {
            beta
          }
        } & type == "BoF"),
        vdm %>% dplyr::filter(beta != {
          {
            beta
          }
        } & type == "GoF")
      )
    gd$index <- toupper(gd$index)
    vdc$index <- toupper(vdc$index)
    vdm$index <- toupper(vdm$index)
    gd <- gd %>% tidyr::drop_na()
    gp <- suppressWarnings(
      ggplot2::ggplot(gd) +
        ggplot2::geom_histogram(
          ggplot2::aes(x = fit, fill = Model),
          alpha = .3,
          position = "identity",
          bins = 30
        ) +
        ggplot2::labs(
          title = "Histograms of simulated fit indices",
          subtitle = "Correct and misspecificed models",
          caption = stringr::str_wrap(
            paste0(
              "Vertical lines illustrate the ",
              alpha,
              " quantiles for correct models (GoF)
                          or ",
              1 - alpha,
              " quantiles for correct models (BoF) in blue and the ",
              1 - beta,
              " quantiles for misspecfiied models (GoF) or ",
              beta,
              " quantiles for misspecified models (BoF) in orange"
            ), width = 40),
          x = "Fit",
          y = "Frequency"
        ) +
        ggplot2::facet_wrap(index ~ ., ncol = 2, scales = "free") +
        ggplot2::geom_vline(ggplot2::aes(xintercept = qc), vdc, color = "darkblue") +
        ggplot2::geom_vline(ggplot2::aes(xintercept = qm), vdm, color = "darkorange") +
        ggplot2::scale_fill_manual(values = c(
          "correct" = "darkblue",
          "misspecified" = "orange"
        )) +
        ggplot2::theme_linedraw() +
        ggplot2::theme(legend.position = "top")
    )
    return(gp)
  }
