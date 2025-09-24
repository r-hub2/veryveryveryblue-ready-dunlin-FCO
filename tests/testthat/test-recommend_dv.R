#context("recommend")
library(FCO)
library(testthat)

mod <- "
F1 =~ Q5 + Q7 + Q8
F2 =~ Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"

fits.single <- gen_fit(mod1 = mod, x = bb1992, rep = 10)

fits.dv.con <- gen_fit(
  mod1 = mod,
  x = bb1992,
  rep = 10,
  dv = TRUE,
  dv.factors = c("F4", "F5"),
  dv.cutoff = .9
)

fits.dv.merge <- gen_fit(
  mod1 = mod,
  x = bb1992,
  rep = 10,
  dv = TRUE,
  dv.factors = c("F4", "F5"),
  merge.mod = TRUE
)

test_that("Is an error produced when fits is from one model?",
          {
            expect_error(suppressWarnings(recommend_dv(fits = fits.single)))
          })

test_that("Are multiple fit indices supported?",
          {
            out <- suppressWarnings(recommend_dv(
              fits = fits.dv.merge,
              index = c("CFI", "SRMR", "GFI")
            ))
            expect_equal(ncol(out$decisions), 3)
          })

test_that("Is the proper dv technique used? merging?",
          {
            out <- suppressWarnings(recommend_dv(
              fits = fits.dv.merge,
            ))
            expect_equal(out$fit.values$model[2], "merged")
          })

test_that("Is the proper dv technique used? constraining?",
          {
            out <- suppressWarnings(recommend_dv(
              fits = fits.dv.con,
            ))
            expect_equal(out$fit.values$model[2], "constrained")
          })

test_that("Is the result rounded properly?",
          {
            out <- suppressWarnings(recommend_dv(
              fits = fits.dv.con,
              digits = 5
            ))
            mm <- rep(NA, nrow(out$cutoffs))
            for (i in seq_along(out$cutoffs[,1])) {
              mm[i] <- match(TRUE, round(out$cutoffs[i, 1], 1:5) == out$cutoffs[i, 1])
            }
            expect_equal(max(mm), 5)
          })
