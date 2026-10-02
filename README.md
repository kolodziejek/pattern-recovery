# Pattern recovery

R code for the tables and figures in Sections 4.1-4.3 of
[From Graphical Lasso to Atomic Norms: High-Dimensional Pattern Recovery](https://arxiv.org/abs/2506.13353).

## Run

Requires R 4.0 or later and the `gmp` package. Install it once:

```r
install.packages("gmp")
```

In R or RStudio, set the working directory to the folder containing `run.R`, then run:

```r
source("run.R")
reproduce("tables")
reproduce("figures")
```

Tables are printed and saved as text, separately for **p = 16** and **p = 64**.
Figures are saved as PNG files. Settings are in `config.R`:
**p = 25**, **1,000 replications** per graph/penalty, and perturbation radii
uniform on **[0, 4 delta]**. Covariance deviations use the entrywise maximum
norm for GLASSO and the projection-adapted norm for ell-infinity and tuned
SLOPE. See the paper for the mathematical definitions.

Results are saved in `results_Q/`. Keep the RDS checkpoints to resume an
interrupted run or redraw completed figures without refitting:

```r
reproduce("replot")
```

After changing numerical code or settings, choose a new output folder:

```r
reproduce("figures", output = "results_run02")
```

## License and citation

MIT License; see `LICENSE`. Please cite the accompanying paper (`CITATION.cff`).
