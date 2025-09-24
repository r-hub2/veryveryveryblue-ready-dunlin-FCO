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

test_that("Is the constraint set properly?",
          {
            out <- FCO:::constr_mod(mod, dv.factors = c("F4", "F5"), dv.cutoff = .9)
            out <- lavaanify(out)
            expect_equal(out[which(out$lhs == "F4" & out$rhs == "F5"), "ustart"], .9)
          })

test_that("Is it also set properly when the order is inversed?",
          {
            out <- FCO:::constr_mod(mod, dv.factors = c("F5", "F4"), dv.cutoff = .9)
            out <- lavaanify(out)
            expect_equal(out[which(out$lhs == "F4" & out$rhs == "F5"), "ustart"], .9)
          })

