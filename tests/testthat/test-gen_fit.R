#context("gen_fit")
library(FCO)
library(testthat)
library(lavaan)
library(semTools)

mod <- "
F1 =~ Q5 + Q7 + Q8
F2 =~ Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"

pop.mod <- pop_mod(mod = mod, x = bb1992)$pop.mod

test_that("Is the result fco of gen_fit a list of length rep?",
          {
            out <- gen_fit(mod1 = mod, x = bb1992, rep = 10)
            expect_type(out$fco, "list")
            expect_equal(length(out$fco), 10)
          })

test_that("Is mod1 a model?",
          {
            out <- gen_fit(mod1 = mod, x = bb1992, rep = 10)
            expect_s3_class(lavaanify(out$mod1), "data.frame")
          })

test_that("Is x a data.frame?",
          {
            out <- gen_fit(mod1 = mod, x = bb1992, rep = 10)
            expect_s3_class(out$x, "data.frame")
          })

test_that("Is n expected when x is not given?",
          {
            expect_error(gen_fit(mod1 = mod, n = 500, rep = 10))
          })

test_that("Is an error returned when neither x nor n is not given?",
          {
            expect_error(gen_fit(mod1 = mod, rep = 10))
          })

test_that("Is an error returned when mod is not given?",
          {
            expect_error(gen_fit(x = bb1992, rep = 10))
          })

test_that("Is an error returned when pop.mod is not given?",
          {
            expect_error(gen_fit(n = 500, rep = 10))
          })

test_that("Is standard type NM?",
          {
            out <- gen_fit(mod1 = mod, x = bb1992, rep = 10)
            expect_equal(out$type, "NM")
          })

test_that("Is type HB taken over?",
          {
            out <- gen_fit(
              mod1 = mod,
              x = bb1992,
              rep = 10,
              type = "HB"
            )
            expect_equal(out$type, "HB")
          })

test_that("Is type EM taken over?",
          {
            out <- gen_fit(
              mod1 = mod,
              x = bb1992,
              rep = 10,
              type = "EM"
            )
            expect_equal(out$type, "EM")
          })

mod2 <- "
F1 =~ Q5 + Q7 + Q8 + Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"

test_that("Is the result fco of gen_fit a list of matrices with two rows?",
          {
            out <- gen_fit(
              mod1 = mod,
              mod2 = mod2,
              x = bb1992,
              rep = 10
            )
            expect_type(out$fco, "list")
            expect_equal(nrow(do.call("rbind", out$fco)), 20)
          })

test_that("Is a second model accepted?",
          {
            out <- gen_fit(
              mod1 = mod,
              mod2 = mod2,
              x = bb1992,
              rep = 10
            )
            expect_s3_class(lavaanify(out$mod2), "data.frame")
          })

test_that("Is type EM used when DV is TRUE?",
          {
            out <- gen_fit(
              mod1 = mod,
              x = bb1992,
              rep = 10,
              dv = TRUE
            )
            expect_equal(out$type, "EM")
          })

test_that("Is constraining used when DV is TRUE?",
          {
            out <- gen_fit(
              mod1 = mod,
              x = bb1992,
              rep = 10,
              dv = TRUE
            )
            expect_equal(out$merge.mod, FALSE)
          })

test_that("Is dv.cutoff set to .9 when DV is TRUE?",
          {
            out <- gen_fit(
              mod1 = mod,
              x = bb1992,
              rep = 10,
              dv = TRUE
            )
            expect_equal(out$dv.cutoff, .9)
          })

test_that(
  "Is the first and second factor merged when DV is TRUE, merge.mod is TRUE and dv.factors is missing?",
  {
    out <-
      gen_fit(
        mod1 = mod,
        x = bb1992,
        rep = 10,
        dv = TRUE,
        merge.mod = TRUE
      )
    expect_equal(grepl("F2 =~ ", out$mod2), FALSE)
  }
)

test_that("Is merging used when DV is TRUE, merge.mod is TRUE and dv.factors is provided?",
          {
            out <-
              gen_fit(
                mod1 = mod,
                x = bb1992,
                rep = 10,
                dv = TRUE,
                merge.mod = TRUE,
                dv.factors = c("F2", "F3")
              )
            expect_equal(out$merge.mod, TRUE)
          })

test_that("Is dv.cutoff set to .85 taken over?",
          {
            out <-
              gen_fit(
                mod1 = mod,
                x = bb1992,
                rep = 10,
                dv = TRUE,
                dv.cutoff = .85
              )
            expect_equal(grepl("F2~~0.85*F1", out$mod2), FALSE)
          })

test_that("Is a warning produced when standardized is switched off?",
          {
            expect_warning(gen_fit(
              mod1 = mod,
              x = bb1992,
              rep = 10,
              standardized = FALSE
            ))
          })

test_that("Is kurtosis and skewness unequal to 1 when assume.mvn is FALSE?",
          {
            out <-
              gen_fit(
                mod1 = mod,
                x = bb1992,
                rep = 10,
                assume.mvn = FALSE
              )
            expect_gt(mardiaKurtosis(out$x)[1], 1)
            expect_gt(mardiaSkew(out$x)[1], 1)
          })

test_that("Is the number of cores 1 when multi.core is FALSE?",
          {
            out <-
              gen_fit(
                mod1 = mod,
                x = bb1992,
                rep = 10,
                multi.core = FALSE
              )
            expect_equal(out$cores, 1)
          })

test_that("Is the seed set properly?",
          {
            out1 <-
              gen_fit(
                mod1 = mod,
                x = bb1992,
                rep = 10
              )
            out2 <-
              gen_fit(
                mod1 = mod,
                x = bb1992,
                rep = 10,
                seed = 1111
              )
            expect_equal(out1$seed, out2$seed)
          })
