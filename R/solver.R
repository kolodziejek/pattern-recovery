# Symmetric scaled ADMM for tr(S K) - log det(K) + 2*lambda*J_w(off(K)).
# The diagonal is unpenalized. Matching lower/upper triangle indices preserve
# symmetry; no rounding, pattern projection, or positive-definite repair is used.
solver_control <- function() list(
  rho = 1, max_iterations = 20000L, absolute_tolerance = 1e-13,
  relative_tolerance = 1e-12, stationarity_tolerance = 1e-12,
  gap_tolerance = 1e-10, kkt_tolerance = 1e-6, check_every = 10L)

sorted_l1 <- function(x, w) sum(w * sort(abs(x), decreasing = TRUE))
sorted_l1_dual <- function(x, w) {
  max(cumsum(sort(abs(x), decreasing = TRUE)) / cumsum(w))
}

# Isotonic least squares is the sorted-L1 proximal map after sorting/sign removal.
# stats::isoreg is a compiled routine shipped with R; no grpSLOPE dependency.
prox_sorted_l1 <- function(x, thresholds) {
  assert(length(x) == length(thresholds) && all(is.finite(x)) &&
         all(is.finite(thresholds)) && all(thresholds >= 0) &&
         all(diff(thresholds) <= 0), "Invalid sorted-L1 proximal arguments.")
  if (!length(x) || max(thresholds) == 0) return(x)
  order_abs <- order(abs(x), decreasing = TRUE, method = "radix")
  z <- abs(x[order_abs]) - thresholds
  # Project onto z_1 >= ... >= z_m >= 0.
  fitted <- if (length(z) == 1L) {
    max(z, 0)
  } else {
    pmax(-stats::isoreg(seq_along(z), -z)$yf, 0)
  }
  ans <- numeric(length(x))
  ans[order_abs] <- sign(x[order_abs]) * fitted
  ans
}
prox_penalty_matrix <- function(A, thresholds, layout = sym_layout(nrow(A))) {
  validate_symmetric(A, "proximal input", tolerance = 1e-10)
  out <- A
  v <- prox_sorted_l1(A[layout$lower_index], thresholds)
  out[layout$lower_index] <- v
  out[layout$upper_index] <- v  # Explicit matching transpose indices.
  out
}
precision_objective <- function(K, S, lambda, w, layout) {
  ld <- logdet_spd(K)
  if (!is.finite(ld)) return(Inf)
  sum(S * K) - ld + 2*lambda*sorted_l1(K[layout$lower_index], w)
}

# Feasible dual: W = S + lambda*Pi, diag(Pi)=0, J_w^*(off(Pi)) <= 1.
# Dual objective p + log det(W). A duality gap is computed, not inferred from
# changes between iterates. No data or iterate is projected to positive definite.
solution_diagnostics <- function(K, S, lambda, w, layout) {
  p <- nrow(K)
  bad <- list(valid = FALSE, objective = Inf, dual_objective = NA_real_,
              gap = Inf, relative_gap = Inf, dual_violation = Inf,
              contact_relative = Inf, stationarity = Inf,
              diagonal_residual = Inf, min_eigenvalue = NA_real_)
  CK <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(CK)) return(bad)
  inverse <- chol2inv(CK)
  T <- symmetrize(inverse - S)
  x <- K[layout$lower_index]
  z <- T[layout$lower_index] / lambda
  nd <- sorted_l1_dual(z, w)
  jx <- sorted_l1(x, w)
  contact <- abs(jx - sum(x*z)) / max(jx, 1e-14)
  station <- max(max_abs(diag(T)),
    max_abs(x - prox_sorted_l1(x + T[layout$lower_index], lambda*w))) /
      max(1, max_abs(K), max_abs(S))
  Pi <- lift_off(z / (max(1, nd) * (1 + 1e-12)), p)
  scale <- 1
  for (attempt in 0:60) {
    W <- symmetrize(S + lambda*scale*Pi)
    CW <- tryCatch(chol(W), error = function(e) NULL)
    if (!is.null(CW)) break
    scale <- scale/2
  }
  if (is.null(CW)) return(bad)
  primal <- sum(S*K) - 2*sum(log(diag(CK))) + 2*lambda*jx
  dual <- p + 2*sum(log(diag(CW)))
  gap <- primal - dual
  # A materially negative gap is a failure, not something to clip away.
  if (gap < -1e-10 * max(1, abs(primal), abs(dual))) return(bad)
  gap <- max(0, gap)
  list(valid = TRUE, objective = primal, dual_objective = dual,
       gap = gap, relative_gap = gap/max(1, abs(primal), abs(dual)),
       dual_violation = max(0, nd - 1), contact_relative = contact,
       stationarity = station,
       diagonal_residual = max_abs(diag(T)),
       min_eigenvalue = min_eigen(K))
}

fit_precision <- function(S, lambda, w, control = solver_control()) {
  defaults <- solver_control()
  assert(is.list(control) && length(control) == length(defaults) &&
         !anyDuplicated(names(control)) && setequal(names(control), names(defaults)),
         "Supply the complete solver_control() list; unknown controls are not supported.")
  assert(identical(as.numeric(control$rho), 1), "This reproduction uses fixed ADMM rho = 1.")
  positive_fields <- c("rho", "max_iterations", "absolute_tolerance",
    "relative_tolerance", "kkt_tolerance", "stationarity_tolerance",
    "gap_tolerance", "check_every")
  assert(all(vapply(control[positive_fields], function(x)
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0, logical(1))),
    "Solver tolerances, iteration limits, check frequency, and rho must be positive.")
  assert(control$max_iterations == as.integer(control$max_iterations) &&
         control$check_every == as.integer(control$check_every),
         "Iteration limits and check frequency must be integers.")
  validate_symmetric(S, "observed covariance")
  spd_chol(S, "observed covariance")
  p <- nrow(S); layout <- sym_layout(p)
  validate_weights(w, layout$m)
  assert(length(lambda) == 1L && is.finite(lambda) && lambda >= 0,
         "lambda must be finite and nonnegative.")
  if (lambda == 0) {
    K <- spd_inverse(S)
    return(list(precision = K, X = K, converged = TRUE, iterations = 0L,
      primal_residual = 0, dual_residual = 0, epsilon_primal = 0,
      epsilon_dual = 0, rho_final = control$rho,
      diagnostics = list(valid = TRUE, objective = p + logdet_spd(S),
       dual_objective = p + logdet_spd(S), gap = 0, relative_gap = 0,
       dual_violation = 0, contact_relative = 0, stationarity = 0,
       diagonal_residual = 0, min_eigenvalue = min_eigen(K))))
  }
  rho <- control$rho
  assert(is.finite(rho) && rho > 0, "ADMM rho must be positive.")
  # Entirely data-based initialization, without K_star or pattern information.
  Y <- spd_inverse(S)
  U <- matrix(0, p, p)
  converged <- FALSE
  diag_result <- NULL
  for (iteration in seq_len(control$max_iterations)) {
    A <- symmetrize(rho * (Y - U) - S)
    e <- eigen(A, symmetric = TRUE)
    root <- sqrt(e$values^2 + 4*rho)
    values <- numeric(p)
    nonnegative <- e$values >= 0
    values[nonnegative] <- (e$values[nonnegative] + root[nonnegative])/(2*rho)
    # Avoid cancellation for large negative eigenvalues.
    values[!nonnegative] <- 2/(root[!nonnegative] - e$values[!nonnegative])
    X <- tcrossprod(sweep(e$vectors, 2L, sqrt(values), "*"))
    Yold <- Y
    Y <- prox_penalty_matrix(symmetrize(X + U), lambda*w/rho, layout)
    U <- symmetrize(U + X - Y)  # SCALED dual update: no extra factor rho.
    primal <- frobenius(X - Y)
    dual <- rho * frobenius(Y - Yold)
    ep <- p * control$absolute_tolerance + control$relative_tolerance *
      max(frobenius(X), frobenius(Y))
    ed <- p * control$absolute_tolerance + control$relative_tolerance *
      frobenius(rho * U)
    residuals_ok <- primal <= ep && dual <= ed
    check <- iteration %% control$check_every == 0L || iteration == control$max_iterations
    if (residuals_ok && check) {
      diag_result <- solution_diagnostics(Y, S, lambda, w, layout)
      converged <- isTRUE(diag_result$valid) &&
        diag_result$relative_gap <= control$gap_tolerance &&
        diag_result$dual_violation <= control$kkt_tolerance &&
        diag_result$contact_relative <= control$kkt_tolerance &&
        diag_result$stationarity <= control$stationarity_tolerance
      if (converged) break
    }

  }
  # Recompute diagnostics at the returned iterate even after a failed stop.
  diag_result <- solution_diagnostics(Y, S, lambda, w, layout)
  list(precision = Y, X = X, converged = converged,
       iterations = iteration, primal_residual = primal, dual_residual = dual,
       epsilon_primal = ep, epsilon_dual = ed, rho_final = rho,
       diagnostics = diag_result)
}

# Compare with the population's target face. GLASSO tests signed support;
# tuned SLOPE also tests its active cluster; ell-infinity tests its maximal
# signed cluster without incorrectly requiring the outside entries to vanish.
pattern_recovered <- function(K, spec, tolerance) {
  assert(is.numeric(tolerance) && length(tolerance) == 1L &&
    is.finite(tolerance) && tolerance > 0, "Pattern tolerance must be positive.")
  x <- K[spec$layout$lower_index]
  if (any(!is.finite(x))) return(FALSE)
  active <- x[spec$active]*spec$signs
  outside <- abs(x[spec$inactive])
  sign_ok <- all(active > tolerance)
  if (spec$pattern == "sign") return(sign_ok && max_abs(outside) <= tolerance)
  tied <- max(active)-min(active) <= tolerance
  if (spec$pattern == "sparse_cluster")
    return(sign_ok && tied && max_abs(outside) <= tolerance)
  assert(spec$pattern == "max_cluster", "Unknown target pattern.")
  separated <- !length(outside) || max(outside) < min(active)-tolerance
  sign_ok && tied && separated
}

# One-sided numerical rejection of a tiny ell-infinity cluster split. The
# subgradient is supported ONLY on exact maximizers of the stored estimate.
# This does not certify recovery or alter K. Margin 10 is a floating-point
# robustness policy, not an interval-arithmetic theorem certificate.
ellinf_split_evidence <- function(K, S, lambda, spec) {
  assert(spec$model == "ellinf" && spec$pattern == "max_cluster" &&
    identical(as.numeric(spec$w), c(1, rep(0, spec$m-1L))),
    "Split evidence requires the unscaled ell-infinity penalty.")
  spd_chol(K, "estimated precision"); spd_chol(S, "observed covariance")
  assert(identical(dim(K),dim(S)) && all(K == t(K)) &&
    length(lambda) == 1L && is.finite(lambda) && lambda > 0,
    "Split evidence requires symmetric matching matrices and positive lambda.")
  x <- K[spec$layout$lower_index]
  active <- x[spec$active]*spec$signs
  spread <- max(active)-min(active)
  distance <- sqrt(2*sum((active-mean(active))^2))
  out <- list(supported_failure=FALSE, residual_F=NA_real_, bound_F=NA_real_,
    active_spread=spread, distance_F=distance, numerical_margin=10,
    reason="one_sided_test_inconclusive")
  if (!length(x) || max(abs(x)) == 0) return(out)
  J <- which(abs(x) == max(abs(x)))
  inverse <- spd_inverse(K)
  z <- sign(x[J])*(inverse-S)[spec$layout$lower_index[J]]/lambda
  # Euclidean simplex projection, shifted first to avoid cancellation from a
  # common offset; its nonnegative masses sum to one on the exact support J.
  z <- z-max(z)
  ordered <- sort(z,decreasing=TRUE); sums <- cumsum(ordered)
  rho <- max(which(ordered-(sums-1)/seq_along(ordered) > 0))
  masses <- pmax(z-(sums[rho]-1)/rho,0)
  masses <- masses/sum(masses)
  masses[which.max(masses)] <- masses[which.max(masses)] + 1-sum(masses)
  assert(all(masses >= 0) && abs(sum(masses)-1) <= 64*.Machine$double.eps,
    "Invalid ell-infinity subgradient masses.")
  pi <- numeric(length(x)); pi[J] <- sign(x[J])*masses
  g <- S-inverse+lambda*lift_off(pi,nrow(K))
  gamma <- frobenius(g)
  kappa <- max(eigen(K,symmetric=TRUE,only.values=TRUE)$values)
  out$residual_F <- gamma
  if (!is.finite(kappa*gamma) || kappa*gamma >= 1) {
    out$reason <- "residual_bound_unavailable"
    return(out)
  }
  out$bound_F <- kappa^2*gamma/(1-kappa*gamma)
  out$supported_failure <- spread > 20*out$bound_F || distance > 10*out$bound_F
  if (out$supported_failure) out$reason <- "target_tie_excluded_with_margin10"
  out
}

# The public classification has two categories only. Recovery requires every
# tested tolerance to pass and no accepted one-sided evidence against the tie.
classify_fit <- function(fit, spec, tolerances = c(1e-10,1e-9,1e-8),
                         S = NULL, lambda = NULL) {
  assert(is.numeric(tolerances) && length(tolerances) > 0L &&
    all(is.finite(tolerances)) && all(tolerances > 0),
    "Pattern tolerances must be positive finite values.")
  flags <- rep(NA, length(tolerances)); evidence <- NULL
  converged <- isTRUE(fit$converged)
  if (converged) {
    flags <- vapply(tolerances,function(tol)
      pattern_recovered(fit$precision,spec,tol),logical(1))
    zero_penalty <- is.numeric(lambda) && length(lambda) == 1L &&
      is.finite(lambda) && lambda == 0
    if (spec$model == "ellinf" && !zero_penalty) {
      assert(!is.null(S) && !is.null(lambda),
        "Ell-infinity classification requires the fitted covariance and lambda.")
      evidence <- ellinf_split_evidence(fit$precision,S,lambda,spec)
    }
  }
  recovery <- converged && all(flags) && !isTRUE(evidence$supported_failure)
  list(status=if(recovery) "recovered" else "failure", recovery=recovery,
    flags=flags, converged=converged, residual_evidence=evidence)
}
