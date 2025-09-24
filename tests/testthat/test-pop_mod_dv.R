#context("pop_mod_dv")
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

pop.mod <- pop_mod(mod = mod, x = bb1992, type = "EM")

test_that("Is pop.mod a model?",
          {
            out <- FCO:::pop_mod_dv(pop.mod = pop.mod$pop.mod, dv.factors = c("F4", "F5"), dv.cutoff = .9)
            expect_true(is.data.frame(lavaanify(out)))
          })

test_that("Does it modify the correlation between specified factors?",
          {
            out <- FCO:::pop_mod_dv(pop.mod = pop.mod$pop.mod, dv.factors = c("F4", "F5"), dv.cutoff = .9)
            out <- lavaanify(out)
            expect_equal(out[which(out$lhs == "F4" & out$rhs == "F5"), "ustart"], .9)
          })

test_that("Does it take the first and second factor if not specified otherwise?",
          {
            out <- FCO:::pop_mod_dv(pop.mod = pop.mod$pop.mod, dv.cutoff = .9)
            out <- lavaanify(out)
            expect_equal(out[which(out$lhs == "F1" & out$rhs == "F2"), "ustart"], .9)
          })

test_that("Does it work if the cutoff is identical to a correlation between other factors?",
          {
            out <- FCO:::pop_mod_dv(pop.mod = pop.mod$pop.mod, dv.factors = c("F4", "F5"), dv.cutoff = .625)
            out <- lavaanify(out)
            expect_equal(out[which(out$lhs == "F4" & out$rhs == "F5"), "ustart"], .625)
          })

test_that("Does it work if the the factor order is reversed?",
          {
            out <- FCO:::pop_mod_dv(pop.mod = pop.mod$pop.mod, dv.factors = c("F5", "F4"), dv.cutoff = .9)
            out <- lavaanify(out)
            expect_equal(out[which(out$lhs == "F4" & out$rhs == "F5"), "ustart"], .9)
          })

test_that("Does it give an error if the factor is not in the model?",
          {
            expect_error(FCO:::pop_mod_dv(pop.mod = pop.mod$pop.mod, dv.factors = c("F4", "F9"), dv.cutoff = .9))
          })
