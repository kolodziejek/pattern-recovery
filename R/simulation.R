# Each replication has its own deterministic seed, so resuming changes no draws.
replication_seed <- function(master, p, graph, model, i) {
  h <- as.double(master %% 2147483646)
  for (a in utf8ToInt(paste(master, "main", p, graph, model, i, sep = ":")))
    h <- (h * 131 + a) %% 2147483646
  as.integer(h + 1)
}

# Draw an independent scalar radius uniformly on [0, multiplier*delta],
# then normalize the symmetric random direction IN THE SELECTED NORM.
# This is not a uniform distribution on the norm ball. For ellinf/SLOPE the
# norm is ||.||_Q; GLASSO retains the entrywise maximum norm.
draw_covariance <- function(spec, Q, metric, upper, seed) {
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion",
           sample.kind = "Rejection")
  H <- matrix(0, spec$p, spec$p)
  H[upper.tri(H, diag = TRUE)] <- runif(spec$p * (spec$p + 1L)/2, -1, 1)
  H[lower.tri(H)] <- t(H)[lower.tri(H)]
  magnitude <- runif(1, 0, upper)
  size <- deviation_size(H, Q, metric)
  assert(is.finite(size) && size > 0, "Degenerate covariance direction.")
  spec$Sigma + magnitude * H/size
}

save_checkpoint <- function(result, path) {
  temporary <- paste0(path, ".tmp")
  saveRDS(result, temporary, version = 3)
  assert(file.rename(temporary, path), paste("Cannot save", path))
}

run_case <- function(analysis, graph, cfg, directory, source_hash) {
  spec <- analysis$spec; constants <- analysis$constants; Q <- analysis$Q
  metric <- constants$deviation_norm
  id <- sprintf("p%d_%s_%s", spec$p, graph, spec$model)
  path <- file.path(directory, "raw", paste0(id, ".rds"))
  meta <- list(version = "projection-norm-v2", graph = graph,
    model = spec$model, p = spec$p, N = cfg$n,
    deviation_norm = metric, deviation_description = deviation_name(metric),
    multiplier = cfg$multiplier, seed = cfg$seed,
    pattern_tolerances = cfg$pattern_tolerances, solver = cfg$solver,
    theory = cfg$theory, source_hash = source_hash)
  assert(constants$valid, paste("No applicable recovery certificate for", id))
  message(id, ": covariance deviation = ", deviation_name(metric),
          "; precision error = entrywise maximum norm.")
  # Domain certificate: ||E||_2 <= p*||E||_max <= p*C*N(E),
  # C=1 for linf and C=2 for the projection-adapted norm.
  # No draws are rejected, truncated, or repaired.
  upper <- cfg$multiplier * constants$delta
  pd_lower_bound <- min_eigen(spec$Sigma) - spec$p * constants$norm_to_max * upper
  assert(pd_lower_bound > 0, paste("Sampling domain is not certified positive definite:", id))
  if (file.exists(path)) {
    result <- readRDS(path)
    assert(identical(result$meta, meta) && identical(result$constants, constants),
      paste("Settings or numerical code changed. Use a new output directory:", id))
  } else {
    zero <- fit_precision(spec$Sigma, constants$lambda, spec$w, cfg$solver)
    label <- classify_fit(zero, spec, cfg$pattern_tolerances,
                          S = spec$Sigma, lambda = constants$lambda)
    assert(zero$converged && label$recovery, paste("Population recovery check failed:", id))
    result <- list(meta = meta, constants = constants, population = spec,
      projection = Q, pd_lower_bound = pd_lower_bound, zero_fit = zero,
      records = data.frame(), complete = FALSE)
    save_checkpoint(result, path)
  }
  if (isTRUE(result$complete)) return(result)
  rows <- if (nrow(result$records)) split(result$records, seq_len(nrow(result$records))) else list()
  done <- length(rows)
  if (done < cfg$n) for (i in seq.int(done + 1L, cfg$n)) {
    seed <- replication_seed(cfg$seed, spec$p, graph, spec$model, i)
    S <- draw_covariance(spec, Q, metric, upper, seed)
    fit <- fit_precision(S, constants$lambda, spec$w, cfg$solver)
    if (!fit$converged) {
      result$records <- if (length(rows)) do.call(rbind, rows) else data.frame()
      save_checkpoint(result, path)
      saveRDS(list(replication = i, seed = seed, covariance = S, fit = fit),
              file.path(directory, "raw", paste0(id, "_failed.rds")))
      stop("Numerical fit not accepted: ", id, " replication ", i,
           ". Diagnostics saved; no recovery label assigned.", call. = FALSE)
    }
    label <- classify_fit(fit, spec, cfg$pattern_tolerances, S, constants$lambda)
    d <- fit$diagnostics
    rows[[i]] <- data.frame(replication = i, seed = seed,
      deviation_norm = metric,
      covariance_error = deviation_size(S - spec$Sigma, Q, metric),
      covariance_max_error = max_abs(S - spec$Sigma),
      covariance_min_eigenvalue = min_eigen(S),
      precision_error = max_abs(fit$precision - spec$K),
      recovered = label$recovery, converged = fit$converged,
      iterations = fit$iterations, primal_residual = fit$primal_residual,
      dual_residual = fit$dual_residual, relative_gap = d$relative_gap,
      stationarity = d$stationarity, dual_violation = d$dual_violation,
      contact_relative = d$contact_relative,
      split_excluded = isTRUE(label$residual_evidence$supported_failure),
      tolerance_sensitive = length(unique(label$flags)) > 1L)
    for (j in seq_along(label$flags)) rows[[i]][[paste0("recovered_tol_", j)]] <- label$flags[j]
    if (i %% cfg$checkpoint_every == 0L || i == cfg$n) {
      result$records <- do.call(rbind, rows)
      save_checkpoint(result, path)
      message(sprintf("%s: %d/%d", id, i, cfg$n))
    }
  }
  result$complete <- TRUE
  save_checkpoint(result, path)
  result
}

run_simulations <- function(cfg, directory, root) {
  assert(cfg$n >= 1 && cfg$n == as.integer(cfg$n) && cfg$multiplier > 0,
         "n must be a positive integer and multiplier must be positive.")
  dir.create(file.path(directory, "raw"), recursive = TRUE, showWarnings = FALSE)
  files <- c("config.R", "R/matrix.R", "R/geometry.R", "R/deviation.R",
             "R/solver.R", "R/simulation.R")
  hashes <- setNames(unname(tools::md5sum(file.path(root, files))), files)
  # Leave an interactive caller's random-number state unchanged.
  old_kind <- RNGkind()
  old_seed <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  results <- list()
  for (graph in cfg$graphs) {
    message("Computing constants: ", graph, ", p = ", cfg$p)
    a <- analyze_graph(cfg$p, graph, cfg$models, cfg$theory, cfg$deviation_norm)
    for (model in cfg$models) results[[paste(graph, model, sep = "_")]] <-
      run_case(a$models[[model]], graph, cfg, directory, hashes)
  }
  results
}

read_results <- function(directory) {
  paths <- list.files(file.path(directory, "raw"), pattern = "^p[0-9]+_.*_(glasso|ellinf|slope)\\.rds$",
                      full.names = TRUE)
  assert(length(paths) == 12L, "Expected twelve saved graph/penalty cases; run the figures command first.")
  results <- lapply(paths, readRDS)
  assert(all(vapply(results, function(r)
    identical(r$meta$version, "projection-norm-v2"), logical(1))),
    "These are older-format results. Generate new simulations in results_Q first.")
  assert(all(vapply(results, function(r) isTRUE(r$complete) && nrow(r$records) == r$meta$N,
                    logical(1))), "The saved simulations are not complete.")
  results
}
