#' Helper function that constraints a given model for discriminant validity testing (constraining)

#' @param mod A lavaan model to be constrained.
#' @param dv.factors Names of the factors to be considered. Must be equal to 2. If missing (the default), the first and second factor of the model are selected.
#' @param dv.cutoff Critical correlation assumed to be a cutoff for discriminant validity testing.
#' For example, based on Rönkkö & Cho (2020), a cutoff of .9 indicates a severe issue in discriminant validity between the selected factors. Cutoffs between .8 and 1 are recommended.
#' The function returns a warning, if the cutoff is below .8.
#' @keywords internal
#' @return An alternative model with the cutoff as correlation between the selected factors. This model can be used in function gen_fit to generate flexible cutoffs. It will be automatically generated in gen_fit, if not provided manually.
#' @noRd
constr_mod <- function(mod, dv.factors = NULL, dv.cutoff) {
  mod <- gsub(" ", "", mod)
  pms <- unlist(strsplit(mod, "\n"))
  vars <- gsub("=~.*", "", pms[grep("=~", pms, fixed = TRUE)])
  if (!is.null(dv.factors))
    mf <- vars[match(dv.factors, vars)]
  if (is.null(dv.factors))
    mf <- vars[c(1, 2)]
  if (any(is.na(mf)))
    stop("At least one of the factors to be merged is not a factor in the model. Please revise.")
  lhs <- grep(paste0(mf[2], "~~"), pms)
  rhs <- grep(paste0("*" , mf[1]), pms)
  mt <- match(lhs, rhs)
  if (length(mt) == 0) {
    lhs2 <- grep(paste0(mf[1], "~~"), pms)
    rhs2 <- grep(paste0("*" , mf[2]), pms)
    mt <- match(lhs2, rhs2)
  }
  if (length(mt) == 0) {
    mod2 <-
      paste(mod, paste0(mf[2], "~~", dv.cutoff, "*", mf[1]), sep = "\n")
  }
  if (length(mt) != 0) {
    if (exists("lhs2")) {
      pms[stats::na.omit(c(lhs2[mt], rhs2[mt]))] <-
        paste0(mf[2], "~~", dv.cutoff, "*", mf[1])
    }
    else {
      pms[stats::na.omit(c(lhs[mt], rhs[mt]))] <-
        paste0(mf[2], "~~", dv.cutoff, "*", mf[1])
    }
    mod2 <- paste(pms, sep = "\n", collapse = "\n")
  }
  return(mod2)
}
