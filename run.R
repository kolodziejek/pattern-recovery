# Terminal: Rscript run.R [tables|figures|replot] [output-directory]
# R/RStudio: source("run.R"); reproduce("tables"); reproduce("figures")
# Requires the gmp package for exact-zero certification.
.source_file <- if (sys.nframe() > 0L) sys.frame(1L)$ofile else
  sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1L])
.project_root <- dirname(normalizePath(.source_file))
for (.file in c("R/matrix.R", "R/geometry.R", "R/deviation.R", "R/solver.R", "config.R",
               "R/simulation.R", "R/reporting.R")) source(file.path(.project_root, .file))

reproduce <- function(task = "tables", output = file.path(.project_root, "results_Q"),
                      cfg = reproduction_config()) {
  task <- match.arg(task, c("tables", "figures", "replot"))
  assert(all(cfg$models %in% names(cfg$deviation_norm)), "Missing deviation-norm settings.")
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  if (task == "tables") {
    analyses <- list()
    for (p in cfg$table_dimensions) for (graph in cfg$graphs) {
      message("Table constants: ", graph, ", p = ", p)
      analyses[[paste(p, graph, sep = "_")]] <- analyze_graph(
        p, graph, cfg$models, cfg$theory, cfg$deviation_norm)
    }
    saveRDS(analyses, file.path(output, "table_constants.rds"), version = 3)
    write_tables(analyses, output)
  }
  if (task %in% c("figures", "replot")) {
    results <- if (task == "replot") read_results(output) else
      run_simulations(cfg, output, .project_root)
    summary <- plot_results(results, output)
    write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
    write_experiment_notes(results, summary, output)
    violations <- sum(vapply(results, function(r)
      sum(!r$records$recovered & r$records$covariance_error < r$constants$delta), integer(1)))
    if (violations) warning(violations, " numerical non-recoveries strictly below delta; inspect raw diagnostics.")
  }
  writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
  message("Output: ", normalizePath(output))
  invisible(output)
}

if (sys.nframe() == 0L) {
  .args <- commandArgs(trailingOnly = TRUE)
  if (length(.args) > 1L) reproduce(.args[1L], .args[2L]) else
    reproduce(if (length(.args)) .args[1L] else "tables")
}
