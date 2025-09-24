#context("pop_mod")
library(FCO)
library(testthat)

mod <- "
F1 =~ Q5 + Q7 + Q8
F2 =~ Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"

test_that("Is pop.mod a model?",
          {
            out <- pop_mod(mod = mod, x = bb1992)
            expect_true(is.data.frame(lavaanify(out$pop.mod)))
          })

test_that("Is mod a model?",
          {
            out <- pop_mod(mod = mod, x = bb1992)
            expect_true(is.data.frame(lavaanify(out$mod)))
          })

test_that("Is the standard type NM?",
          {
            out <- pop_mod(mod = mod, x = bb1992)
            expect_equal(out$type, "NM")
          })

test_that("Is type HB taken over?",
          {
            out <- pop_mod(mod = mod, x = bb1992, type = "HB")
            expect_equal(out$type, "HB")
          })

test_that("Is type EM taken over?",
          {
            out <- pop_mod(mod = mod, x = bb1992, type = "EM")
            expect_equal(out$type, "EM")
          })

test_that("Is a warning produced when standardized is switched off?",
          {
            expect_warning(pop_mod(
              mod = mod,
              x = bb1992,
              standardized = FALSE
            ))
          })

test_that("Is an alternative setting for afl and aco taken over?",
          {
            out <- pop_mod(mod = mod, x = bb1992, afl = .8, aco = .5)
            expect_equal(out$type, "NM")
            expect_equal(grepl("0.8", out$pop.mod), TRUE)
            expect_equal(grepl("0.5", out$pop.mod), TRUE)
          })

test_that("Are afl and aco overwritten when the type is not NM?",
          {
            out <- pop_mod(mod = mod, x = bb1992, type = "EM", afl = .666, aco = .555)
            expect_equal(out$type, "EM")
            expect_equal(grepl("0.666", out$pop.mod), FALSE)
            expect_equal(grepl("0.555", out$pop.mod), FALSE)
          })
