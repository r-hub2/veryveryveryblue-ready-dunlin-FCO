#' Obtain factor loadings as in Hu & Bentler (1999)
#'
#' @param x Vector of length (1:10) containing the loadings per factor
#' @keywords internal
#' @return Vector of factor loadings for popluation model
#' @noRd
hb_load <- function(x) {
  checkmate::assertInteger(
    length(x),
    lower = 1,
    upper = 10,
    any.missing = FALSE,
    null.ok = FALSE
  )
  if (length(x) == 1)
    y <- 1
  else {
    l <- c(.7, .7, .7, .7, .75, .75, .8, .8, .8, .8)
    m <-
      list(
        c(5:6),
        c(4, 6, 8),
        c(4:7),
        c(3:4, 6, 8:9),
        c(3:6, 8:9),
        c(2:4, 6, 8:10),
        c(2:6, 8:10),
        c(1:4, 6, 7:10),
        c(1:10)
      )
    y <- l[m[[length(x) - 1]]]
  }
  return(y)
}
