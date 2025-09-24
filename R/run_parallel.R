#' Helper function for parallel wrapping
#' @param fun see parallel
#' @param seeds see gen_fit2
#' @param cores see gen_fit2
#' @param method see parallel
#' @keywords internal
#' @importFrom foreach %dopar%
#' @noRd
run_parallel <- function(fun, seeds, cores, method = "foreach", ...) {
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  res <- suppressMessages({
    foreach::foreach(i = seq_along(seeds), .combine = rbind) %dopar% {
      fun(seed = seeds[i], ...)
    }
  })
  res
}
