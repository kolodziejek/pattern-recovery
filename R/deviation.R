# Auxiliary norm for covariance deviations. The penalty itself never changes.
# Q is the canonical population projection stored in orthonormal svec coordinates.
# These routines evaluate symmetric matrices; Q_norms() separately evaluates the
# full ambient operator norms needed for the theorem.
deviation_name <- function(metric) switch(metric,
  linf = "entrywise maximum norm ||.||_infinity",
  Q = "projection-adapted maximum norm ||.||_Q",
  stop("Unknown deviation norm: ", metric, call. = FALSE))

deviation_size <- function(E, Q, metric) {
  if (metric == "linf") return(max_abs(E))
  assert(metric == "Q", "Deviation norm must be 'linf' or 'Q'.")
  v <- svec(E, Q$layout)
  qv <- apply_Q(Q, v)
  max(max_abs(qv / Q$layout$scale),
      max_abs((v - qv) / Q$layout$scale))
}

# Exact restricted comparison c_Q for the two experimental cluster patterns.
# SLOPE: restricted dual unit ball has endpoints +/- mean(w[1:k])*sign_vec.
# ell-infinity: its vertices are +/- sign_vec/k and +/- e_j, j inactive.
# Maximize the norm on these finitely many vertices. No random approximation.
projection_comparison_constant <- function(spec, Q, block_size = 128L) {
  c <- spec$face_multiplier * deviation_size(spec$shadow, Q, "Q")
  if (spec$space == "sparse_cluster") return(c)
  assert(spec$space == "max_cluster",
         "The Q-norm specialization supports ell-infinity and tuned SLOPE clusters.")
  indices <- spec$p + spec$inactive
  if (length(indices)) for (first in seq.int(1L, length(indices), by = block_size)) {
    cols <- indices[first:min(length(indices), first + block_size - 1L)]
    A <- if (ncol(Q$L)) Q$L %*% t(Q$R[cols, , drop = FALSE]) else
      matrix(0, Q$layout$n, length(cols))
    locations <- cbind(cols, seq_along(cols))
    if (Q$identity) A[locations] <- A[locations] + 1
    # svec(D e_j) = sqrt(2)*e_(p+j); divide rowwise to undo svec scaling.
    A <- sweep(sqrt(2) * A, 1L, Q$layout$scale, "/")
    remainder <- -A
    remainder[locations] <- remainder[locations] + 1
    c <- max(c, max_abs(A), max_abs(remainder))
  }
  c
}

# For n_Q, q=b=1 exactly. Face directions and the irrepresentability vector
# lie in ker(Q), and Gamma*vec(Delta) lies in range(Q) for Delta in M_star.
# Hence tau, alpha, and the existing zeta LOWER bound are unchanged.
# The residual constant is NOT 1 and is NOT copied unchanged from linf:
# ||v||_infinity <= 2*n_Q(v) <= 2*L*||v||_infinity,
# L=max(||Q||_infinity,||I-Q||_infinity), gives the valid choice
# eta_Q = 2*L*eta_linf (not a claim that this is the sharp constant).
deviation_constants <- function(spec, Q, norms, metric, block_size = 128L) {
  assert(metric %in% c("linf", "Q"), "Unknown auxiliary deviation norm.")
  L <- max(norms$Qnorm, norms$IQnorm)
  if (metric == "linf") return(list(
    cdiamond = spec$cdiamond, eta = spec$eta,
    Qnorm = norms$Qnorm, IQnorm = norms$IQnorm,
    norm_to_max = 1, eta_description = "valid entrywise-maximum bound"))
  list(cdiamond = projection_comparison_constant(spec, Q, block_size),
    eta = 2 * L * spec$eta, Qnorm = 1, IQnorm = 1, norm_to_max = 2,
    eta_description = "valid upper bound: 2*max(Qnorm_linf,IQnorm_linf)*eta_linf")
}
