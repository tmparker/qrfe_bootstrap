# QRFE bootstrap

Replication files for "Bootstrap inference for panel data quantile regression"

Uploaded 2022-01-10

## Contents

This repository contains files that were used to conduct the Monte Carlo
simulation experiments in the paper.

- The main file used for simulation is `sim.R`.
- The files `run_simulations.sh`, `glue.R`, `missing.R` and `make_tables.R` are
  specific to the distributed computing environment that was used (a
  SHARCNET/Compute Canada system).  They were used to serially farm `sim.R` out
  to many CPUs, "glue" the parts back together, and then extract the information
  from the glued files.
- `se_experiment.R` extracted information about the performance of variance
  estimates computed using the bootstrap from the data files.
- `data/` is a directory containing the results of the simulations for five
  different data generating processes considered in the paper.
- `tables/` is a directory containing the Latex-format tables that appear in the
  paper.

## Basic usage

If you are using R, then the following pseudo-code will compute one bootstrap
coefficient vector.  It assumes you have response vector `y`, `n*T` by `p`
covariate matrix `x` and `n*T` by `n` matrix of unit indicators `ind`, assuming
the panel is balanced (`n` observations with `T` time-series observations per
unit).  Suppose you would like to estimate a quantile regression at quantile
level `Tau`.  The commands below create `nboot` different bootstrap coefficient
estimates for the `p` covariates and you can calculate whatever else you wish
from there.

```
library(quantreg)
w <- rexp(n * nboot, 1)
U <- matrix(rep(w, each = T), n * T, nboot)
bcoef <- boot.rq.wxy(cbind(x, ind), y, U, tau = Tau)[, 1:p]
```
