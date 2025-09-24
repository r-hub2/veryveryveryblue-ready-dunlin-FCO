#context("flex_co")
library(FCO)
library(testthat)

mod <- "
F1 =~ Q5 + Q7 + Q8
F2 =~ Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"
mod2 <- "
F1 =~ Q5 + Q7 + Q8 + Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"

fits <- gen_fit(mod1 = mod, x = bb1992, rep = 10)
fits2 <- gen_fit(mod1 = mod, mod2 = mod2, x = bb1992, rep = 10)

test_that("Is the function warning because of low replications?",
          {
            expect_warning(flex_co(
              fits = fits,
              index = c("CFI", "SRMR", "GFI"),
              alpha.lev = .05,
              gof = c(T, F, T)
            ))
          })

test_that("Is the function accepting multiple indices?",
          {
            out <- suppressWarnings(flex_co(
              fits = fits,
              index = c("CFI", "SRMR", "GFI"),
              alpha.lev = .05,
              gof = c(T, F, T)
            ))
            expect_equal(length(out$cutoff), 3)
            expect_equal(length(out$index), 3)
            expect_equal(length(out$gof), 3)
          })

test_that("Is the function reversing quantile definition based on gof?",
          {
            out <- suppressWarnings(flex_co(
              fits = fits,
              index = c("CFI", "SRMR", "GFI"),
              alpha.lev = .05,
              gof = c(F, T, F)
            ))
            expect_gt(out$cutoff[1], .97)
            expect_lt(out$cutoff[2], .35)
            expect_gt(out$cutoff[3], .96)
          })

test_that("Is the function rejecting multiple levels of uncertainty?",
          {
            expect_error(suppressWarnings(flex_co(
              fits = fits,
              index = c("CFI", "SRMR", "GFI"),
              alpha.lev = c(.05, .10),
              gof = c(T, F, T)
            )))
          })

test_that("Is the function returning values for two models?",
          {
            out <- suppressWarnings(flex_co(
              fits = fits2,
              index = c("CFI", "SRMR", "GFI"),
              alpha.lev = .05,
              gof = c(T, F, T)
            ))
            expect_equal(nrow(out$cutoff), 2)
            expect_equal(length(out$difference), 3)
          })
