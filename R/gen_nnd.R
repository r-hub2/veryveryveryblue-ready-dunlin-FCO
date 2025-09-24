#' Helper function to obtain intermediate matrix for mixed data.types
#' @param x see gen_fit2
#' @param data.types ...
#' @param seed ...
#' @keywords internal
#' @noRd
gen_nnd <- function(x, data.types, seed = seed) {
  if (ncol(x) != length(data.types))
    stop(
      "Please check the data.types vector. It should be of the same length as the number of manifest variables in your model."
    )
  no.count <- length(which(data.types == "C"))
  no.binary <- length(which(data.types == "B"))
  no.ordinal <- length(which(data.types == "O"))
  no.normal <- length(which(data.types == "N"))
  tx <- x[, c(
    which(data.types == "C"),
    which(data.types == "B"),
    which(data.types == "O"),
    which(data.types == "N")
  )]
  lam.count <- apply(x[, which(data.types == "C")], 2, median, na.rm = T)
  p.binary <- apply(x[, which(data.types == "B")], 2, mean, na.rm = T)
  p.ordinal <- unname(apply(x[, which(data.types == "O")], 2, function(x) {
    prop.table(table(x))
  }))
  p.ordinal <- lapply(p.ordinal, cumsum)
  for (i in 1:length(p.ordinal)) {
    p.ordinal[[i]] <- p.ordinal[[i]][p.ordinal[[i]] != 0 &
                                       p.ordinal[[i]] != 1]
  }
  mean.normal <- apply(x[, which(data.types == "N")], 2, mean, na.rm = T)
  var.normal <- apply(x[, which(data.types == "N")], 2, var, na.rm = T)
  rx <- cor(tx, method = "spearman")
  set.seed(seed)
  intm <- PoisBinOrdNor::intermat(
    no_pois = no.count,
    no_bin = no.binary,
    no_ord = no.ordinal,
    no_norm = no.normal,
    corr_mat = rx,
    lam_vec = lam.count,
    prop_vec_bin = p.binary,
    prop_vec_ord = p.ordinal,
    nor_mean = mean.normal,
    nor_var = var.normal
  )
  cx <- PoisBinOrdNor::genPBONdata(
    n = nrow(x),
    no_pois = no.count,
    no_bin = no.binary,
    no_ord = no.ordinal,
    no_norm = no.normal,
    inter.mat = intm,
    lamvec = lam.count,
    prop_vec_bin = p.binary,
    prop_vec_ord = p.ordinal,
    nor.mean = mean.normal,
    nor.var = var.normal
  )
  cx <- data.frame(cx$data)
  names(cx) <- names(tx)
  cx <- as_tibble(cx)
  cx <- cx[, names(x)]
  return(cx)
}
