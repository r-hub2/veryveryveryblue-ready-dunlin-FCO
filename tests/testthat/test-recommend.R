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
mod2 <- "
F1 =~ Q5 + Q7 + Q8 + Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"

fits <- gen_fit(mod1 = mod, x = bb1992, rep = 10)
fits2 <- gen_fit(
  mod1 = mod,
  mod2 = mod2,
  x = bb1992,
  rep = 10
)

test_that("Are the settings for purpose and focus following the recommendations by Mai et al.? #1",
          {
            out <- suppressWarnings(recommend(
              fits = fits,
              purpose = "novel",
              focus = "CFA",
              override = FALSE,
            ))
            expect_equal(rownames(out$recommended), "SRMR")
          })

test_that("Are the settings for purpose and focus following the recommendations by Mai et al.? #2",
          {
            out <- suppressWarnings(recommend(
              fits = fits,
              purpose = "established",
              focus = "CFA",
              override = FALSE,
            ))
            expect_equal(rownames(out$recommended), "CFI")
          })

test_that("Are the settings for purpose and focus following the recommendations by Mai et al.? #3",
          {
            out <- suppressWarnings(recommend(
              fits = fits,
              purpose = "novel",
              focus = "structural",
              override = FALSE,
            ))
            expect_equal(rownames(out$recommended), "SRMR")
          })

test_that("Are the settings for purpose and focus following the recommendations by Mai et al.? #4",
          {
            out <- suppressWarnings(
              recommend(
                fits = fits,
                purpose = "established",
                focus = "structural",
                override = FALSE,
              )
            )
            expect_equal(rownames(out$recommended), "SRMR")
          })

test_that("Is an error produced when override is chosen without index?",
          {
            expect_error(suppressWarnings(recommend(fits = fits,
                                                    override = T, )))
          })

test_that("Is an error produced when fits is from two models?",
          {
            expect_error(suppressWarnings(recommend(fits = fits2)))
          })

test_that("Is the result rounded properly?",
          {
            out <- suppressWarnings(
              recommend(
                fits = fits,
                purpose = "established",
                focus = "structural",
                override = FALSE,
                digits = 5
              )
            )
            mm <- rep(NA, nrow(out$cutoffs))
            for (i in seq_along(out$cutoffs[,1])) {
              mm[i] <- match(TRUE, round(out$cutoffs[i, 1], 1:5) == out$cutoffs[i, 1])
            }
            expect_equal(max(mm), 5)
          })

test_that("Are the indices properly used when override is TRUE?",
          {
            out <- suppressWarnings(
              recommend(
                fits = fits,
                purpose = "established",
                focus = "structural",
                override = TRUE,
                index = c("GFI", "TLI")
              )
            )
            expect_equal(names(out$cutoffs), c("GFI", "TLI"))
          })
