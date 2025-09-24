#context("check_std")
library(FCO)
library(lavaan)
library(testthat)

mod <- "
F1 =~ Q5 + Q7 + Q8
F2 =~ Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"

fit <- cfa(mod, bb1992, auto.fix.first = F, std.lv = T)

test_that("Is .9 really .9 in a standardized model?",
          {
            out <- FCO:::check_std(fit, dv.factors = c("F4", "F5"), dv.cutoff = .9)
            expect_equal(out, .9)
          })
