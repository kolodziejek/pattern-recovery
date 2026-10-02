# Deterministic constants for the one-nonzero-magnitude graph families.
# Full vec uses column-major order. Projectors are the canonical orthogonal
# projectors onto a subspace of Sym(p), not projectors onto nonsymmetric arrays.

validate_weight_dimension <- function(m) {
  assert(is.numeric(m) && length(m) == 1L && is.finite(m) &&
         m >= 1 && m == floor(m), "m must be a positive integer scalar.")
  invisible(TRUE)
}
validate_active_count <- function(k, m) {
  validate_weight_dimension(m)
  assert(is.numeric(k) && length(k) == 1L && is.finite(k) &&
         k >= 1 && k <= m && k == floor(k),
         "k must be an integer scalar between 1 and the weight-vector length m.")
  invisible(TRUE)
}
validate_weights <- function(w, m) {
  validate_weight_dimension(m)
  assert(is.numeric(w) && is.null(dim(w)) && length(w) == m && all(is.finite(w)),
         "Weights must be a finite numeric vector of the required length m (p*(p-1)/2 for matrices).")
  assert(all(w >= 0) && w[1] > 0 && all(diff(w) <= 0),
         "Weights must be nonincreasing, nonnegative, with w[1] > 0.")
  invisible(TRUE)
}
make_weights <- function(model, m, k) {
  model <- match.arg(model, c("glasso", "ellinf", "slope"))
  validate_active_count(k, m)
  if (model == "glasso") w <- rep(1, m)
  if (model == "ellinf") w <- c(1, rep(0, m - 1L))
  if (model == "slope") {
    if (k == 1L) {
      w <- rep(1, m)
    } else if (k == m) {
      j <- k %/% 2L
      w <- c(rep(1, j), if (k %% 2L) 0.5, rep(0, j))
    } else {
      j <- k %/% 2L
      w <- c(rep(1, j), if (k %% 2L) 2/3, rep(1/3, m - j - k %% 2L))
    }
  }
  validate_weights(w, m)
  w
}

# Exact infinity-ball radius for the specialized sparse one-cluster subspace.
# A present active block must be nonconstant when k >= 2, and a present
# tail must contain a positive weight. Degenerate blocks change the pattern
# subspace and must use an explicit geometry branch, not just a different tau.
# For the active permutahedron: every subset-size inequality is checked.
# This computes the actual face radius, including repeated weights.
cluster_face_radius <- function(w, k) {
  m <- length(w)
  validate_weights(w, m)
  validate_active_count(k, m)
  a <- w[seq_len(k)]
  if (k == 1L) {
    active <- Inf
  } else {
    assert(max(a) > min(a),
           "A constant active weight block has a different pattern subspace.")
    j <- seq_len(k - 1L)
    active <- min((cumsum(a)[j] - j * mean(a)) / pmin(j, k - j))
  }
  if (k < m) assert(w[k + 1L] > 0,
    paste("An all-zero inactive weight block has a different pattern subspace;",
          "use the ell-infinity branch for (1,0,...,0), not the sparse-cluster projector."))
  tail <- if (k < m) mean(w[(k + 1L):m]) else Inf
  min(active, tail)
}

population_spec <- function(K, model, tolerance = 1e-12) {
  model <- match.arg(model, c("glasso", "ellinf", "slope"))
  assert(is.numeric(tolerance) && length(tolerance) == 1L &&
         is.finite(tolerance) && tolerance >= 0,
         "Structure tolerance must be a finite nonnegative scalar.")
  spd_chol(K, "K_star")
  p <- nrow(K); lay <- sym_layout(p)
  x <- K[lay$lower_index]
  active <- which(abs(x) > tolerance)
  inactive <- setdiff(seq_len(lay$m), active)
  assert(length(active) > 0L, "The all-zero off-diagonal pattern is outside this experiment suite.")
  rho <- mean(abs(x[active]))
  assert(max_abs(abs(x[active]) - rho) <= tolerance * max(1, rho),
         "These geometry routines require one common nonzero absolute edge value.")
  assert(max_abs(x[inactive]) <= tolerance, "Inconsistent support tolerance.")
  k <- length(active); m <- lay$m
  w <- make_weights(model, m, k)
  signs <- sign(x[active])
  sign_vec <- numeric(m); sign_vec[active] <- signs
  if (model == "glasso" || (model == "slope" && k == 1L)) {
    space <- "support"; pattern <- "sign"
    tau <- if (k == m) Inf else 1
    face_multiplier <- 1; cdiamond <- 1
    distance <- rho
  } else if (model == "ellinf") {
    space <- "max_cluster"; pattern <- "max_cluster"
    tau <- if (k == 1L) Inf else 1 / k
    face_multiplier <- 1/k; cdiamond <- 1
    distance <- if (k < m) rho/2 else rho
  } else {
    space <- "sparse_cluster"; pattern <- "sparse_cluster"
    tau <- cluster_face_radius(w, k)
    face_multiplier <- mean(w[seq_len(k)])
    cdiamond <- face_multiplier
    distance <- rho
  }
  shadow <- lift_off(sign_vec, p)
  d <- max(rowSums(abs(shadow))) + 1
  Sigma <- spd_inverse(K)
  eta_sparse <- min(d * row_norm(Sigma) * row_norm(K)^2, sum(abs(K)))
  eta <- if (space == "max_cluster") sum(abs(K)) else eta_sparse
  list(K = K, Sigma = Sigma, p = p, m = m, k = k, layout = lay,
       active = active, inactive = inactive, signs = signs,
       sign_vec = sign_vec, shadow = shadow, rho = rho, d = d,
       w = w, space = space, pattern = pattern, tau = tau,
       face_multiplier = face_multiplier, cdiamond = cdiamond,
       eta = eta, zeta_lower = distance / row_norm(K)^2, model = model)
}

unit_columns <- function(n, indices) {
  out <- matrix(0, n, length(indices))
  if (length(indices)) out[cbind(indices, seq_along(indices))] <- 1
  out
}
active_contrasts <- function(spec) {
  n <- spec$layout$n; k <- spec$k
  out <- matrix(0, n, max(0L, k - 1L))
  if (k > 1L) for (j in seq_len(k - 1L)) {
    v <- c(rep(1, j), -j) / sqrt(j * (j + 1))
    rows <- spec$p + spec$active[seq_len(j + 1L)]
    out[rows, j] <- spec$signs[seq_len(j + 1L)] * v
  }
  out
}
pattern_basis <- function(spec, complement = FALSE) {
  p <- spec$p; n <- spec$layout$n
  if (spec$space == "support") {
    indices <- if (complement) p + spec$inactive else c(seq_len(p), p + spec$active)
    return(unit_columns(n, indices))
  }
  if (complement) {
    contrasts <- active_contrasts(spec)
    if (spec$space == "max_cluster") return(contrasts)
    return(cbind(unit_columns(n, p + spec$inactive), contrasts))
  }
  cluster <- numeric(n)
  cluster[p + spec$active] <- spec$signs / sqrt(spec$k)
  free <- if (spec$space == "max_cluster") c(seq_len(p), p + spec$inactive) else seq_len(p)
  cbind(unit_columns(n, free), cluster)
}
pattern_rank <- function(spec) switch(spec$space,
  support = spec$p + spec$k,
  sparse_cluster = spec$p + 1L,
  max_cluster = spec$p + length(spec$inactive) + 1L)

hessian_times_basis <- function(A, B, layout) {
  out <- matrix(0, nrow(B), ncol(B))
  if (ncol(B)) for (j in seq_len(ncol(B))) {
    Y <- A %*% smat(B[, j], layout) %*% A
    out[, j] <- svec(symmetrize(Y), layout)
  }
  out
}

# Store Q_s as identity + L R', or simply L R'.
# H is the log-det Hessian on orthonormal symmetric coordinates.
# Direct: Q_s = H B (B'H B)^(-1) B'.
# Complement: Q_s = I - C (C'H^(-1)C)^(-1) C'H^(-1).
# Use whichever requires the smaller Cholesky factorization.
make_Q <- function(spec, control) {
  n <- spec$layout$n
  rank <- pattern_rank(spec)
  complement <- (n - rank < rank)
  B <- pattern_basis(spec, complement)
  if (ncol(B) == 0L) return(list(identity = TRUE,
    L = matrix(0, n, 0), R = matrix(0, n, 0), rank = rank,
    factor_dimension = 0L, rcond = 1, layout = spec$layout))
  A <- if (complement) spec$K else spec$Sigma
  HB <- hessian_times_basis(A, B, spec$layout)
  small <- symmetrize(crossprod(B, HB))
  rc <- rcond(small)
  assert(is.finite(rc) && rc >= control$condition_floor,
         sprintf("Ill-conditioned restricted Hessian (rcond %.3g); no pseudoinverse fallback.", rc))
  inv_small <- chol2inv(spd_chol(small, "restricted Hessian"))
  W <- HB %*% inv_small
  list(identity = complement, L = if (complement) -B else W,
       R = if (complement) W else B, rank = rank,
       factor_dimension = ncol(B), rcond = rc, layout = spec$layout)
}
apply_Q <- function(q, x) {
  y <- if (ncol(q$L)) drop(q$L %*% crossprod(q$R, x)) else numeric(length(x))
  if (q$identity) y <- y + x
  y
}
Q_block <- function(q, rows) {
  out <- if (ncol(q$L)) {
    q$L[rows, , drop = FALSE] %*% t(q$R)
  } else {
    matrix(0, length(rows), q$layout$n)
  }
  if (q$identity) out[cbind(seq_along(rows), rows)] <-
      out[cbind(seq_along(rows), rows)] + 1
  out
}

# Exact induced infinity norms of the p^2-dimensional canonical lift.
# In particular ||I-Q|| includes skew-symmetric inputs, as the manuscript's
# norm is defined on all R^(p^2), not just on symmetric perturbations.
Q_norms <- function(q, block_size = 128L) {
  n <- q$layout$n; p <- q$layout$p; scale <- q$layout$scale
  qnorm <- iqnorm <- iq_symmetric <- 0
  for (start in seq.int(1L, n, by = block_size)) {
    rows <- start:min(n, start + block_size - 1L)
    A <- Q_block(q, rows)
    dd <- A[cbind(seq_along(rows), rows)]
    A <- sweep(sweep(A, 2L, scale, "*"), 1L, scale[rows], "/")
    total <- rowSums(abs(A))
    rest <- pmax(0, total - abs(dd))
    diag_contribution <- abs(1 - dd)
    off <- rows > p
    diag_contribution[off] <- abs(1 - dd[off]/2) + abs(dd[off]/2)
    qnorm <- max(qnorm, total)
    iqnorm <- max(iqnorm, rest + diag_contribution)
    iq_symmetric <- max(iq_symmetric, rest + abs(1 - dd))
  }
  list(Qnorm = qnorm, IQnorm = iqnorm, IQnorm_symmetric_domain = iq_symmetric)
}

recovery_constants <- function(spec, control, Q = NULL, norms = NULL,
                               deviation = "linf") {
  if (is.null(Q)) Q <- make_Q(spec, control)
  if (is.null(norms)) norms <- Q_norms(Q, control$block_size)
  selected <- deviation_constants(spec, Q, norms, deviation, control$block_size)
  cc <- selected$cdiamond; eta <- selected$eta
  v <- svec(spec$shadow, spec$layout)
  off_sign <- (apply_Q(Q, v) - v) / spec$layout$scale
  unscaled_lhs <- if (irrep_is_exact_zero(spec)) {
    0
  } else {
    max_abs(off_sign)
  }
  lhs <- spec$face_multiplier * unscaled_lhs
  alpha_max <- if (is.infinite(spec$tau)) 1 else 1 - lhs / spec$tau
  alpha <- if (is.infinite(spec$tau)) NA_real_ else
    max(0, alpha_max) * (1 - control$alpha_safety)
  good <- is.infinite(spec$tau) || (is.finite(alpha) && alpha > 0 && spec$tau > 0)
  out <- c(list(model = spec$model, p = spec$p, k = spec$k, d = spec$d,
    space = spec$space, deviation_norm = deviation,
    deviation_description = deviation_name(deviation),
    tau = spec$tau, cdiamond = cc, eta = eta, zeta_lower = spec$zeta_lower,
    cdiamond_linf = spec$cdiamond, eta_linf = spec$eta,
    Qnorm_linf = norms$Qnorm, IQnorm_linf = norms$IQnorm,
    norm_to_max = selected$norm_to_max, eta_description = selected$eta_description,
    face_multiplier = spec$face_multiplier,
    unscaled_lhs = unscaled_lhs, irrep_lhs = lhs,
    irrep_ratio = if (is.finite(spec$tau)) lhs / spec$tau else 0,
    alpha_max = alpha_max, alpha = alpha,
    valid = good, status = if (good) "valid" else "irrepresentability_not_satisfied",
    factor_dimension = Q$factor_dimension, restricted_rcond = Q$rcond),
    list(Qnorm = selected$Qnorm, IQnorm = selected$IQnorm))
  out$M <- out$r_unconstrained <- out$r <- out$delta <- out$lambda <- NA_real_
  out$error_bound <- out$transformed_error_bound <- NA_real_
  out$radius_truncated <- NA
  if (!good) return(out)
  qn <- selected$Qnorm; iqn <- selected$IQnorm
  singleton <- is.infinite(spec$tau)
  a <- if (singleton) NA_real_ else alpha * spec$tau
  M <- if (singleton) 1/qn else a / (qn * (iqn * cc + a))
  # Stable evaluation when M is tiny.
  root <- sqrt(1 + M)
  r0 <- M / (root * (root + 1))
  radius_cap <- (1 - control$radius_safety) * eta * spec$zeta_lower
  r <- min(r0, radius_cap)
  budget <- (M * r - r*r/(1 - r)) / eta
  assert(is.finite(budget) && budget > 0, "The recovery budget is not strictly positive.")
  if (singleton) {
    lambda <- control$singleton_penalty_fraction * budget / cc
    delta <- budget - lambda * cc
  } else {
    lambda <- r * iqn / (eta * qn * (iqn * cc + a))
    delta <- budget
  }
  assert(lambda > 0 && is.finite(lambda) && r < 1, "Invalid tuning/radius.")
  # Check both sufficient inequalities, rather than trusting the closed form.
  primal_budget <- r/(eta*qn) - lambda*cc - r*r/(eta*(1-r))
  dual_budget <- if (singleton) Inf else lambda*a/iqn - r*r/(eta*(1-r))
  assert(delta <= min(primal_budget, dual_budget) + 1e-11*max(delta, 1e-12),
         "Closed-form threshold failed its defining inequality check.")
  out$M <- M; out$r_unconstrained <- r0; out$r <- r
  out$delta <- delta; out$lambda <- lambda
  out$radius_truncated <- r < r0
  out$transformed_error_bound <- r / eta
  out$error_bound <- row_norm(spec$K)^2 * r / eta
  out
}

# Classical constants use ORDERED support, exactly as in Gamma_{S,S}.
# This is separate from the symmetric-space computation above.
graph_components <- function(mask) {
  p <- nrow(mask); unseen <- seq_len(p); out <- list()
  while (length(unseen)) {
    component <- unseen[1]; frontier <- component
    while (length(frontier)) {
      neighbors <- which(colSums(mask[frontier, , drop = FALSE]) > 0)
      new <- setdiff(neighbors, component)
      component <- union(component, new); frontier <- new
    }
    out[[length(out) + 1L]] <- sort(component)
    unseen <- setdiff(unseen, component)
  }
  out
}
ravikumar_constants <- function(K, control) {
  p <- nrow(K); S <- spd_inverse(K)
  mask <- abs(K) > control$structure_tolerance
  components <- graph_components(mask)
  cliques <- all(vapply(components, function(a) all(mask[a, a, drop = FALSE]), logical(1)))
  if (cliques) {
    # Disconnected complete blocks: the support-restricted inverse Hessian
    # is blockwise K_C (x) K_C and the off-support coupling is exactly zero.
    kg <- max(vapply(components, function(a) row_norm(K[a, a, drop = FALSE])^2, numeric(1)))
    lhs <- 0; rc <- NA_real_
  } else {
    idx <- which(mask); other <- which(!mask)
    ri <- (idx - 1L) %% p + 1L; ci <- (idx - 1L) %/% p + 1L
    H <- S[ri, ri, drop = FALSE] * S[ci, ci, drop = FALSE]
    rc <- rcond(H)
    assert(is.finite(rc) && rc >= control$condition_floor,
           "Classical ordered-support Hessian is ill-conditioned.")
    inv <- chol2inv(spd_chol(H, "ordered-support Hessian"))
    kg <- row_norm(inv); lhs <- 0
    if (length(other)) for (start in seq.int(1L, length(other), by = control$block_size)) {
      block <- other[start:min(length(other), start + control$block_size - 1L)]
      rb <- (block - 1L) %% p + 1L; cb <- (block - 1L) %/% p + 1L
      Hcross <- S[rb, ri, drop = FALSE] * S[cb, ci, drop = FALSE]
      lhs <- max(lhs, row_norm(Hcross %*% inv))
    }
  }
  amax <- 1 - lhs
  a <- max(0, amax) * (1 - control$alpha_safety)
  ks <- row_norm(S); d <- max(rowSums(mask))
  beta <- min(abs(K[lower.tri(K) & mask]))
  delta <- if (a > 0) 1 / (6 * ks^3 * kg^2 * d * (1 + 8/a)^2) else NA_real_
  sign_cap <- if (a > 0) beta / (2*(1+8/a)*kg) else NA_real_
  delta_sign <- if (a > 0) min(delta, (1 - control$radius_safety)*sign_cap) else NA_real_
  list(classical_lhs = lhs, classical_alpha_max = amax, classical_alpha = a,
       kappa_Sigma = ks, kappa_Gamma = kg,
       delta_R_corrected = delta, delta_R_sign = delta_sign,
       classical_beta_min_ok = if (a > 0) delta < sign_cap else FALSE,
       lambda_R_corrected = if (a > 0) 8*delta/a else NA_real_,
       lambda_R_sign = if (a > 0) 8*delta_sign/a else NA_real_,
       classical_rcond = rc)
}

# Deterministic table/experiment analysis, sharing identical pattern projectors.
analyze_graph <- function(p, graph, models, theory_control,
                          deviation_norm = c(glasso = "linf", ellinf = "Q", slope = "Q")) {
  K <- make_precision(graph, p)
  geometries <- list(); out <- list()
  for (model in models) {
    spec <- population_spec(K, model, theory_control$structure_tolerance)
    key <- spec$space
    if (is.null(geometries[[key]])) {
      q <- make_Q(spec, theory_control)
      geometries[[key]] <- list(q = q, norms = Q_norms(q, theory_control$block_size))
    }
    a <- geometries[[key]]
    metric <- unname(deviation_norm[model])
    assert(length(metric) == 1L && !is.na(metric), "Missing deviation norm for this penalty.")
    constants <- recovery_constants(spec, theory_control, a$q, a$norms, metric)
    out[[model]] <- list(spec = spec, constants = constants, Q = a$q)
  }
  list(graph = graph, p = p, models = out,
       classical = ravikumar_constants(K, theory_control))
}


# Certify a zero irrepresentability vector using exact rational arithmetic.
# Since F belongs to the pattern space:
#   (Q - P) vec(F) = 0  iff  K_star F K_star belongs to that space.
# No numerical tolerance or small-value cutoff is used here.
irrep_is_exact_zero <- function(spec) {
  if (!requireNamespace("gmp", quietly = TRUE)) {
    stop("Install the 'gmp' package for exact zero certification.",
         call. = FALSE)
  }
  
  K <- gmp::as.bigq(spec$K)
  F <- gmp::as.bigq(spec$shadow)
  H <- gmp::`%*%`(gmp::`%*%`(K, F), K)
  
  # Extract entries individually to preserve their exact rational values.
  ij <- spec$layout$lower
  off <- lapply(seq_len(nrow(ij)), function(j) {
    H[ij[j, 1L], ij[j, 2L]]
  })
  
  inactive_zero <- all(vapply(spec$inactive, function(j) {
    isTRUE(off[[j]] == 0)
  }, logical(1)))
  
  if (spec$space == "support") {
    return(inactive_zero)
  }
  
  signed_active <- lapply(seq_along(spec$active), function(j) {
    off[[spec$active[j]]] * spec$signs[j]
  })
  
  active_tied <- all(vapply(signed_active, function(x) {
    isTRUE(x == signed_active[[1L]])
  }, logical(1)))
  
  switch(
    spec$space,
    max_cluster = active_tied,
    sparse_cluster = inactive_zero && active_tied,
    stop("Unsupported pattern space.", call. = FALSE)
  )
}
