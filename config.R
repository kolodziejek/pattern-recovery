# Main reproduction settings. No contributed R packages are required.
theory_control <- function() list(
  structure_tolerance = 1e-12, alpha_safety = 1e-10,
  radius_safety = 1e-8, condition_floor = 1e-12,
  block_size = 128L, singleton_penalty_fraction = 0.25)

reproduction_config <- function() list(
  table_dimensions = c(16L, 64L), p = 25L, n = 1000L, multiplier = 4,
  graphs = c("chain", "hub", "grid", "dense"),
  models = c("glasso", "ellinf", "slope"), seed = 260056L,
  # Auxiliary COVARIANCE-deviation norms, not penalty norms:
  # GLASSO: ||vec(E)||_infinity.
  # ell-infinity and tuned SLOPE: ||vec(E)||_Q =
  # max(||Q vec(E)||_infinity, ||(I-Q) vec(E)||_infinity).
  # Each penalty uses its own population pattern projector Q. The selected
  # norm determines delta, lambda, the perturbation magnitude and the x-axis.
  # The y-axis remains max(abs(K_hat-K_star)) in every figure.
  deviation_norm = c(glasso = "linf", ellinf = "Q", slope = "Q"),
  pattern_tolerances = c(1e-10, 1e-9, 1e-8), checkpoint_every = 100L,
  theory = theory_control(), solver = solver_control())
