#' Internal function that takes a model and determines an alternative population model with a fixed correlation between selected factors for discriminant validity testing (constraining)
#'
#' @param pop.mod A population model, potentially from function pop_mod.
#' @param dv.factors (same as constr_mod)
#' @param dv.cutoff (same as constr_mod)
#' @keywords internal
#' @return An alternative population model with the cutoff as correlation between the selected factors. This population model can be used in function gen_fit to generate flexible cutoffs.
#' @noRd
pop_mod_dv <-
  function(pop.mod,
           dv.factors = NULL,
           dv.cutoff) {
    pop.mod <- gsub(" ", "", pop.mod)
    pms <- unlist(strsplit(pop.mod, "\n"))
    vars <- gsub("=~.*", "", pms[grep("=~", pms, fixed = TRUE)])
    if (!is.null(dv.factors))
      mf <- vars[match(dv.factors, vars)]
    if (is.null(dv.factors))
      mf <- vars[c(1, 2)]
    if (any(is.na(mf)))
      stop("At least one of the factors to be constrained is not a factor in the model. Please revise.")
    lhs <- grep(paste0(mf[2], "~~"), pms)
    rhs <- grep(paste0("*" , mf[1]), pms)
    #Maybe the order of factors is inversed in pms?
    if (!any(lhs %in% rhs)) {
      lhs2 <- grep(paste0(mf[1], "~~"), pms)
      rhs2 <- grep(paste0("*" , mf[2]), pms)
      if (any(lhs2 %in% rhs2)) {
        cvr <- pms[rhs2[stats::na.omit(match(lhs2, rhs2))]]
        cvr <-
          as.numeric(unlist(strsplit(gsub(
            ".*~~", "", cvr
          ), "*", fixed = TRUE))[1])
        pop.mod2 <-
          gsub(
            paste0(mf[1], "~~", cvr, "*", mf[2]),
            paste0(mf[1], "~~", dv.cutoff, "*", mf[2]),
            pop.mod,
            fixed = TRUE
          )
      }
    }
    if (any(lhs %in% rhs)) {
      cvr <- pms[rhs[stats::na.omit(match(lhs, rhs))]]
      cvr <-
        as.numeric(unlist(strsplit(gsub(".*~~", "", cvr), "*", fixed = TRUE))[1])
      pop.mod2 <-
        gsub(
          paste0(mf[2], "~~", cvr, "*", mf[1]),
          paste0(mf[2], "~~", dv.cutoff, "*", mf[1]),
          pop.mod,
          fixed = TRUE
        )
    }
    return(pop.mod2)
  }
