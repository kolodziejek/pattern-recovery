# Tables and figures use the computed constants and saved observations directly.

table_number <- function(x) {
  stopifnot(is.numeric(x), all(is.finite(x)))
  ifelse(x == 0, "0", formatC(x, format = "e", digits = 1L))
}

write_tables <- function(analyses, directory) {
  directory <- file.path(directory, "tables")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  graphs <- c("chain", "hub", "grid", "dense")
  labels <- c("Chain graph", "Hub graph", "Grid", "Dense graph")
  tables <- list()
  for (p in c(16L, 64L)) {
    rows <- lapply(graphs, function(graph) {
      a <- Filter(function(a) a$p == p && a$graph == graph, analyses)
      stopifnot(length(a) == 1L)
      a[[1L]]
    })
    irrep <- data.frame("Graph Type" = labels, check.names = FALSE)
    for (model in c("glasso", "ellinf", "slope")) {
      column <- c(glasso = "1", ellinf = "infty", slope = "w")[[model]]
      irrep[[column]] <- table_number(vapply(rows, function(a)
        a$models[[model]]$constants$unscaled_lhs, numeric(1)))
    }
    bounds <- data.frame("Graph Type" = labels,
      delta = table_number(vapply(rows, function(a)
        a$models$glasso$constants$delta, numeric(1))),
      delta_R = table_number(vapply(rows, function(a)
        a$classical$delta_R_corrected, numeric(1))), check.names = FALSE)
    tables[[paste0("irrepresentability_p", p)]] <- irrep
    tables[[paste0("deviation_bounds_p", p)]] <- bounds
  }
  for (name in names(tables)) {
    lines <- c(name, capture.output(print(tables[[name]], row.names = FALSE)))
    cat(paste(lines, collapse = "\n"), "\n\n", sep = "")
    writeLines(lines, file.path(directory, paste0(name, ".txt")))
  }
  norm_note <- c(
    "Auxiliary covariance-deviation norm: GLASSO = ||.||_infinity; ell-infinity and tuned SLOPE = ||.||_Q.",
    "||v||_Q = max(max(abs(Q v)), max(abs(v - Q v))). Q is penalty/population-specific.",
    "The irrepresentability table reports the unscaled face-specific left side, which is unchanged by this norm switch.",
    "The deviation-bound table is a GLASSO comparison; it retains the entrywise maximum norm.")
  # Record the actual configuration, including intentional user overrides.
  choices <- unique(unlist(lapply(analyses, function(a) vapply(a$models, function(b)
    paste(b$spec$model, "=", b$constants$deviation_description), character(1)))))
  norm_note[1L] <- paste("Configured norms:", paste(choices, collapse = "; "))
  cat(paste(norm_note, collapse = "\n"), "\n")
  writeLines(norm_note, file.path(directory, "norms.txt"))
  invisible(tables)
}

# Recovery is decided by the experiment's tolerance and validity checks.
# Plotting must never silently reclassify observations or hide solver failures.
plot_data <- function(result) {
  d <- result$records
  stopifnot(is.data.frame(d), nrow(d) > 0,
    all(c("covariance_error", "precision_error", "recovered", "converged") %in% names(d)))
  if (anyNA(d$converged) || !all(d$converged))
    stop("Cannot publish a figure with unconverged fits; inspect the saved diagnostics.")
  stopifnot(is.logical(d$recovered), !anyNA(d$recovered),
    all(is.finite(d$covariance_error)), all(is.finite(d$precision_error)),
    all(d$covariance_error >= 0), all(d$precision_error >= 0))
  data.frame(x = d$covariance_error, y = d$precision_error, recovered = d$recovered)
}

plot_panel <- function(result, d) {
  delta <- result$constants$delta
  empirical <- if (all(d$recovered)) NA_real_ else min(d$x[!d$recovered])
  subtitle <- if (is.finite(empirical))
    sprintf("Threshold ratio: %.1f", empirical / delta) else "Threshold ratio: not observed"
  # Base R expands the observed range slightly; zero need not be included.
  plot(d$x, d$y, pch = 16, cex = 0.48,
    col = ifelse(d$recovered, "royalblue3", "firebrick3"),
    xlim = range(c(0, d$x, delta)), xlab = "", ylab = "", bty = "l")
  abline(v = delta, col = "grey45", lty = 2, lwd = 1.4)
  if (is.finite(empirical)) abline(v = empirical, col = "darkorange2", lty = 2, lwd = 1.4)
  title(main = graph_title(result$meta$graph, result$meta$p), line = 1.7, cex.main = 1.05)
  mtext(subtitle, side = 3, line = 0.3, cex = 0.9)
  data.frame(graph = result$meta$graph, model = result$meta$model, p = result$meta$p,
    recovery = sum(d$recovered), failure = sum(!d$recovered), delta = delta,
    deviation_norm = result$meta$deviation_norm,
    deviation_description = result$meta$deviation_description,
    empirical_threshold = empirical, threshold_ratio = empirical / delta,
    lambda = result$constants$lambda, multiplier = result$meta$multiplier,
    max_observed_ratio = max(d$x)/delta,
    tolerance_sensitive = sum(result$records$tolerance_sensitive),
    threshold_observed = is.finite(empirical))
}

plot_results <- function(results, directory) {
  directory <- file.path(directory, "figures")
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  graphs <- c("chain", "hub", "grid", "dense")
  summaries <- list()
  for (model in c("glasso", "ellinf", "slope")) {
    panels <- lapply(graphs, function(graph) {
      selected <- Filter(function(r) r$meta$model == model && r$meta$graph == graph, results)
      stopifnot(length(selected) == 1L)
      selected[[1L]]
    })
    data <- lapply(panels, plot_data)
    draw <- function() {
      layout(matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.12))
      par(mar = c(3.1, 4.6, 3.6, 2.2), oma = c(0, 0, 2, 0),
          mgp = c(2, 0.7, 0), las = 1, cex.axis = 0.85)
      summary <- do.call(rbind, Map(plot_panel, panels, data))
      par(mar = rep(0, 4))
      plot.new()
      legend("center", legend = c("Recovery", "Failure", expression(delta), "Observed threshold"),
        col = c("royalblue3", "firebrick3", "grey45", "darkorange2"),
        pch = c(16, 16, NA, NA), lty = c(NA, NA, 2, 2), lwd = c(NA, NA, 1.4, 1.4),
        horiz = TRUE, bty = "n", cex = 0.95)
      mtext(model_title(model), side = 3, outer = TRUE, line = 0.3, font = 2, cex = 1.25)
      summary
    }
    png(file.path(directory, paste0(model, ".png")), width = 1800, height = 1584,
        res = 180)
    summaries[[model]] <- tryCatch(draw(), finally = dev.off())
  }
  invisible(do.call(rbind, summaries))
}

# Norm choice is explicit in the console and saved output, NEVER on the plot.
write_experiment_notes <- function(results, summary, directory) {
  shown <- data.frame(Graph = summary$graph, Penalty = summary$model,
    `Deviation norm` = ifelse(summary$deviation_norm == "Q", "||.||_Q", "||.||_infinity"),
    `Threshold ratio` = ifelse(summary$threshold_observed,
      sprintf("%.1f", summary$threshold_ratio), "not observed"), check.names = FALSE)
  print(shown, row.names = FALSE)
  constants <- do.call(rbind, lapply(results, function(r) data.frame(
    graph = r$meta$graph, model = r$meta$model, p = r$meta$p,
    deviation_norm = r$meta$deviation_norm,
    deviation_description = r$meta$deviation_description,
    tau = r$constants$tau, c = r$constants$cdiamond, eta = r$constants$eta,
    eta_description = r$constants$eta_description,
    Qnorm = r$constants$Qnorm, IQnorm = r$constants$IQnorm,
    Qnorm_linf = r$constants$Qnorm_linf, IQnorm_linf = r$constants$IQnorm_linf,
    alpha = r$constants$alpha, zeta_lower = r$constants$zeta_lower,
    M = r$constants$M, r = r$constants$r, delta = r$constants$delta,
    lambda = r$constants$lambda, pd_lower_bound = r$pd_lower_bound)))
  write.csv(constants, file.path(directory, "simulation_constants.csv"), row.names = FALSE)
  notes <- c("Numerical reproduction settings", "",
    "The covariance-deviation norm and its constants are recorded per case below.",
    "For ||.||_Q: ||v||_Q=max(||Q v||_infinity,||(I-Q)v||_infinity).",
    "Q=Gamma*(P Gamma P)^+ uses the true population and the selected penalty's pattern.",
    "The vertical coordinate is always max(abs(K_hat-K_star)). Neither axis is labelled on the plots.",
    "Threshold ratio = smallest observed numerical non-recovery deviation / delta (rounded to 0.1 on plots).",
    "No observed non-recovery is reported as 'not observed', never as an infinite threshold or a proven lower bound.",
    "Recovery requires agreement at every configured pattern tolerance and no accepted one-sided split evidence.",
    "Tolerance-sensitive observations are numerical non-recoveries under this convention; they are counted separately in summary.csv.",
    "The multiplier is a fixed experimental choice; the range is never adjusted after observing recovery labels.",
    "", capture.output(print(shown, row.names = FALSE)), "",
    vapply(results, function(r) sprintf(
      "%s / %s: %s; p=%d, N=%d, radius range=[0, %.1f*delta], lambda=%.9g, delta=%.9g, PD lower bound=%.9g.",
      r$meta$graph, r$meta$model, r$meta$deviation_description, r$meta$p, r$meta$N,
      r$meta$multiplier, r$constants$lambda, r$constants$delta, r$pd_lower_bound), character(1)))
  writeLines(notes, file.path(directory, "experiment_settings.txt"))
}
