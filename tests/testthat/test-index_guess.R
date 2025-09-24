#context("index_guess")
library(FCO)
library(testthat)

#Index Guess
test_that("Is CFI a GoF?",
          {
            out <- index_guess("CFI")
            expect_equal(out, "GoF")
          })
test_that("Is SRMR a BoF?",
          {
            out <- index_guess("SRMR")
            expect_equal(out, "BoF")
          })
test_that("Is something else?",
          {
            out <- index_guess("Donald")
            expect_equal(out, "not a fit index")
          })
