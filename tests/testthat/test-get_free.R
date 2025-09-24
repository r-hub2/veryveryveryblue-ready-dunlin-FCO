#context("get_free")
library(FCO)
library(lavaan)
library(testthat)

mod <- "
F1 =~ Q5 + Q7 + Q8
F2 =~ Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
F4 ~~ 0.9 * F5
"

fit <- cfa(mod, bb1992, auto.fix.first = F, std.lv = T)
pt <- parTable(fit)


test_that("Is a  constrained correlation kept if constraining is selected?",
          {
            out <- FCO:::get_free(pt, dv.factors = c("F4", "F5"), mode = "constraining")
            out <- lavaanify(out)
            expect_equal(out[which(out$lhs == "F4" & out$rhs == "F5"), "ustart"], .9)
          })

test_that("Is a any constrained correlation omitted if merging is selected?",
          {
            out <- FCO:::get_free(pt, dv.factors = c("F4", "F5"), mode = "merged")
            out <- lavaanify(out)
            expect_equal(nrow(out[which(out$lhs == "F4" & out$rhs == "F5"), ]), 0)
          })


