#' Internal function to get pop_mod for a regression-including model
#' @keywords internal
#' @param mod see gen_fit2
#' @param x ...
#' @param type ...
#' @param es.lam ...
#' @param es.cor ...
#' @param es.f2 ...
#' @noRd
pop_mod_reg <- function(mod, x, type = "EM", es.lam = "low", es.cor = "large", es.f2 = "moderate") {
  #Load propper reference values for NM and HB
  if (type != "EM") {
    steps <- c("low", "moderate", "large")
    if (!es.lam %in% steps) stop("Please use 'low', 'moderate' or 'large' only for es.lam.")
    if (!es.cor %in% steps) stop("Please use 'low', 'moderate' or 'large' only for es.cor.")
    if (!es.f2 %in% steps) stop("Please use 'low', 'moderate' or 'large' only for es.f2.")
    ref.lam <- c(.7, .8, .9)
    names(ref.lam) <- steps
    ref.lam <- ref.lam[which(names(ref.lam) == es.lam)]
    ref.cor <- c(.1, .3, .5)
    names(ref.cor) <- steps
    ref.cor <- ref.cor[which(names(ref.cor) == es.cor)]
    ref.f2 <- c(.02, .15, .35)
    names(ref.f2) <- steps
    ref.f2 <- ref.f2[which(names(ref.f2) == es.f2)]
  }
  if (type == "EM") {
    ref.lam <- NULL; ref.cor <- NULL; ref.f2 <- NULL
  }

  #Model estimation
  fit <- lavaan::sem(
    model = mod,
    data = x,
    std.lv = F,
    auto.fix.first = T
  )
  #Make sure that this is from a model with std.lv = T and auto.fix.first = F
  if (!fit@Options$std.lv | fit@Options$auto.fix.first) {
    fit <- lavaan::sem(model = mod, data = fit@Data, std.lv = T, auto.fix.first = F)
  }
  pt <- lavaan::parTable(fit)

  #Empty objects
  n.ov <- NULL; n.f <- NULL
  lam.mod <- NULL; cov.mod <- NULL; rv.mod <- NULL
  ovv.mod <- NULL; reg.mod <- NULL

  #Observed
  n.ov <-
    unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(rhs))$rhs
  #Factors
  n.f <-
    unique(pt %>% dplyr::filter(op == "=~") %>% dplyr::select(lhs))$lhs

  #Loadings per factor
  lam.list <- lapply(seq_along(n.f), function(i) {
    tl <- pt %>%
      dplyr::filter(op == "=~") %>%
      dplyr::filter(lhs == n.f[i])
    if (nrow(tl) > 1) {
      if (type == "NM") {
        tl$est <- rep(ref.lam, nrow(tl))
      }
      if (type == "HB") {
        if (ref.lam == .7) hbc <- -.05
        if (ref.lam == .8) hbc <- .05
        if (ref.lam == .9) hbc <- .15
        tl$est <- hb_load(tl$est) + hbc
      }
    }
    lam.mod <- paste0(n.f[i], "=~", paste0(round(tl$est, 3), "*", tl$rhs, collapse = "+"))
    list(mod = lam.mod, lam = tl)
  })
  lam.mod <- sapply(lam.list, "[[", "mod")
  lam.2 <- lapply(lam.list, "[[", "lam")

  #Factor covariances aka correlations
  n.cov <- pt %>% dplyr::filter(op == "~~") %>% dplyr::filter(lhs != rhs) %>% dplyr::filter(lhs %in% n.f)
  if (nrow(n.cov > 0)) {
    cov.list <- lapply(seq_along(1:nrow(n.cov)), function(i) {
      tl <- n.cov[i, ]
      if (type == "NM") {
        tl$est <- ref.cor
      }
      if (type == "HB") {
        tl$est <- rep(c(ref.cor - .1, ref.cor, ref.cor + .1), nrow(n.cov))[i]
      }
      cov.mod <- paste0(tl$lhs, "~~", paste0(round(tl$est, 3), "*", tl$rhs))
      list(mod = cov.mod, cov = tl)
    })
    cov.mod <- sapply(cov.list, "[[", "mod")
  }

  #Observed covariances aka residual covariances
  n.rv <- pt %>% dplyr::filter(op == "~~") %>% dplyr::filter(lhs != rhs) %>% dplyr::filter(lhs %in% n.ov)
  if (nrow(n.rv) > 0) {
    rv.list <- lapply(seq_along(1:nrow(n.rv)), function(i) {
      tl <- n.rv[i,]
      if (type == "NM" | type == "HB") {
        tl$est <- 0
      }
      rv.mod <- paste0(tl$lhs, "~~", paste0(round(tl$est, 3), "*", tl$rhs))
      list(mod = rv.mod, rv = tl)
    })
    rv.mod <- sapply(rv.list, "[[", "mod")
  }

  #Factor variances
  n.fv <- pt %>% dplyr::filter(op == "~~") %>% dplyr::filter(lhs == rhs) %>% dplyr::filter(lhs %in% n.f)
  if (nrow(n.fv) > 0) {
    fv.list <- lapply(seq_along(1:nrow(n.fv)), function(i) {
      tl <- n.fv[i,]
      if (type == "NM" | type == "HB") {
        tl$est <- 1
      }
      fv.mod <- paste0(tl$lhs, "~~", paste0(round(tl$est, 3), "*", tl$rhs))
      list(mod = fv.mod, fv = tl)
    })
    fv.mod <- sapply(fv.list, "[[", "mod")
  }

  #Observed variances
  n.ovv <- pt %>% dplyr::filter(op == "~~") %>% dplyr::filter(lhs == rhs) %>% dplyr::filter(lhs %in% n.ov) %>% dplyr::filter(free != 0)
  if (nrow(n.ovv) > 0) {
    tl <- n.ovv
    if (type == "NM" | type == "HB")
      tl$est <- 1 - do.call("rbind", lam.2)$est^2
    ovv.mod <- paste0(tl$lhs, "~~", paste0(round(tl$est, 3), "*", tl$rhs))
  }

  #Regressions
  n.reg <- pt %>% dplyr::filter(op == "~") %>% dplyr::filter(lhs != rhs) %>% dplyr::filter(lhs %in% n.f)
  if (nrow(n.reg) == 0) warning("The model is a CFA. Please use pop_mod, not pop_mod_reg.")
  if (nrow(n.reg) > 0) {
    reg.list <- lapply(seq_along(1:nrow(n.reg)), function(i) {
      tl <- n.reg[i,]
      if (type == "NM" | type == "HB")
        #Cohen from R2 = F2 / (F2 + 1)
        tl$est <- sqrt(ref.f2 / (ref.f2 + 1))
      reg.mod <- paste0(tl$lhs, "~", paste0(round(tl$est, 3), "*", tl$rhs))
      list(mod = reg.mod, reg = tl)
    })
    reg.mod <- sapply(reg.list, "[[", "mod")
  }

  #Factor intercepts
  n.if <- pt %>% dplyr::filter(op == "~1") %>% dplyr::filter(lhs %in% n.f) %>% dplyr::filter(free != 0)
  if (nrow(n.if) > 0) {
    if.list <- lapply(seq_along(1:nrow(n.if)), function(i) {
      tl <- n.if[i,]
      tl$est <- 1
      if.mod <- paste0(tl$lhs, "~1")
      list(mod = if.mod, if.tl = tl)
    })
    if.mod <- sapply(if.list, "[[", "mod")
  }

  #Observed intercepts
  n.iov <- pt %>% dplyr::filter(op == "~1") %>% dplyr::filter(lhs %in% n.ov) %>% dplyr::filter(free != 0)
  if (nrow(n.iov) > 0) {
    iov.list <- lapply(seq_along(1:nrow(n.iov)), function(i) {
      tl <- n.iov[i,]
      tl$est <- 1
      iov.mod <- paste0(tl$lhs, "~1")
      list(mod = iov.mod, iov.tl = tl)
    })
    iov.mod <- sapply(iov.list, "[[", "mod")
  }

  #Resulting population model
  #if and iov omitted so far
  res.mod <- list(lam.mod, cov.mod, rv.mod, fv.mod, ovv.mod, reg.mod) #if.mod, iov.mod
  res.mod <- res.mod[!sapply(res.mod, is.null)]
  pop.mod <- paste0(unlist(res.mod), collapse = "\n")
  res <- list(pop.mod = pop.mod, mod = mod, x = x, type = type, ref.lam = ref.lam, ref.cor = ref.cor, ref.f2 = ref.f2)
  return(res)
}
