# Basic utilities and the four precision-matrix families.
# No contributed R packages are required.
`%||%` <- function(x, y) if (is.null(x)) y else x
max_abs <- function(x) if (length(x)) max(abs(x)) else 0
row_norm <- function(x) if (length(x)) max(rowSums(abs(x))) else 0
frobenius <- function(x) sqrt(sum(x * x))
symmetrize <- function(x) (x + t(x)) / 2
assert <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)

validate_symmetric <- function(x, name = "matrix", tolerance = 1e-12) {
  assert(is.matrix(x) && is.numeric(x) && nrow(x) == ncol(x),
         paste(name, "must be a numeric square matrix."))
  assert(all(is.finite(x)), paste(name, "contains non-finite values."))
  assert(max_abs(x - t(x)) <= tolerance * max(1, max_abs(x)),
         paste(name, "is not symmetric; it will not be silently repaired."))
  invisible(TRUE)
}
spd_chol <- function(x, name = "matrix") {
  validate_symmetric(x, name)
  out <- tryCatch(chol(x), error = function(e) NULL)
  assert(!is.null(out), paste(name, "must be positive definite."))
  out
}
spd_inverse <- function(x, name = "matrix") chol2inv(spd_chol(x, name))
logdet_spd <- function(x) {
  a <- tryCatch(chol(x), error = function(e) NULL)
  if (is.null(a)) return(NA_real_)
  2 * sum(log(diag(a)))
}
min_eigen <- function(x) min(eigen(x, symmetric = TRUE, only.values = TRUE)$values)

make_precision <- function(graph, p) {
  assert(length(p) == 1L && is.finite(p) && p >= 3 && p == as.integer(p),
         "p must be an integer at least 3.")
  p <- as.integer(p)
  graph <- match.arg(graph, c("chain", "hub", "grid", "dense"))
  K <- diag(p)
  if (graph == "chain") {
    diag(K) <- c(25/24, rep(26/24, p - 2L), 25/24)
    K[cbind(seq_len(p - 1L), 2:p)] <- -5/24
    K[cbind(2:p, seq_len(p - 1L))] <- -5/24
  } else if (graph == "hub") {
    q <- 2.5 / (p - 1)
    assert(q < 1, "The manuscript's hub parameter requires p >= 4.")
    diag(K) <- 1 / (1 - q*q)
    K[1, -1] <- K[-1, 1] <- -q / (1 - q*q)
    K[1, 1] <- (1 + (p - 2)*q*q) / (1 - q*q)
  } else if (graph == "grid") {
    side <- as.integer(round(sqrt(p)))
    assert(side * side == p, "The square grid requires p to be a perfect square.")
    for (i in seq_len(side)) for (j in seq_len(side)) {
      a <- (i - 1L) * side + j
      if (i < side) K[a, a + side] <- K[a + side, a] <- 0.1
      if (j < side) K[a, a + 1L] <- K[a + 1L, a] <- 0.1
    }
  } else {
    K[seq_len(p - 1L), seq_len(p - 1L)] <- 0.1
    diag(K) <- 1
  }
  spd_chol(K, paste(graph, "population precision"))
  K
}
graph_title <- function(graph, p) switch(graph,
  chain = "Chain graph", hub = "Hub graph", dense = "Dense graph",
  grid = sprintf("Grid %d x %d", as.integer(sqrt(p)), as.integer(sqrt(p))))
model_title <- function(model) switch(model,
  glasso = "GLASSO", ellinf = "ell-infinity", slope = "Tuned SLOPE")

# Orthonormal coordinates on Sym(p): diagonal, then sqrt(2)*lower triangle.
# These coordinates are NOT the manuscript's vec_+; convert explicitly.
sym_layout <- function(p) {
  lower <- which(lower.tri(matrix(FALSE, p, p)), arr.ind = TRUE)
  list(p = p, m = nrow(lower), n = p + nrow(lower), lower = lower,
       lower_index = lower[, 1] + (lower[, 2] - 1L)*p,
       upper_index = lower[, 2] + (lower[, 1] - 1L)*p,
       scale = c(rep(1, p), rep(sqrt(2), nrow(lower))))
}
svec <- function(x, layout = sym_layout(nrow(x))) {
  c(diag(x), sqrt(2) * x[layout$lower_index])
}
smat <- function(x, layout) {
  ans <- matrix(0, layout$p, layout$p)
  diag(ans) <- x[seq_len(layout$p)]
  v <- x[layout$p + seq_len(layout$m)] / sqrt(2)
  ans[layout$lower_index] <- v
  ans[layout$upper_index] <- v
  ans
}
lift_off <- function(v, p) {
  lay <- sym_layout(p)
  assert(length(v) == lay$m, "Incorrect off-diagonal vector length.")
  a <- matrix(0, p, p)
  a[lay$lower_index] <- v
  a[lay$upper_index] <- v
  a
}
