# This script investigates the performance of two smoothed quantile regression
# estimators compared to the standard estimator (i.e., using the quantile
# regression check function for estimation).
# This script generates data, estimates models and conducts inference over
# simulation repetitions.  Comparisons are processed in another script.

rm(list=ls())
library(quantreg)
library(rlecuyer)
library(conquer)

source("../../common.R")
source("../common_smooth.R")

# 100 instances, each instance runs 10 repetitions.
# Simulation design 1 means arg should run 0-99
arg <- as.numeric(commandArgs(trailingOnly = TRUE))
nstream <- 100 # number of parallel random number streams
dnum <- ceiling((arg + 1) / nstream) # design number
pnum <- arg %% nstream + 1
stream_names <- paste0("st", 1:nstream)
.lec.SetPackageSeed(6:1)
.lec.CreateStream(stream_names) # Seed info stored in list .lec.Random.seed.table
# 200 different per-part RNG seeds set below around the simulation run.

sim_run <- function(n, T, des_num, reps, nboot, taus = 1:3 / 4, alpha = 0.1) {
  ntaus <- length(taus)
  al2 <- alpha / 2
  z_al2 <- qnorm(1 - al2)
  d <- kronecker(diag(n), rep(1, T))
  des <- design(des_num)
  nml <- des$nml
  het <- des$het
  dep <- des$dep
  dyn <- des$dyn
  # Static and dynamic models are generated much differently
  if (dyn) {
    true_beta <- beta_true(taus, be = 0.4) # just 0.4s for dynamic model
  } else {
    true_beta <- beta_true(taus, nml = nml, het = het, dep = dep)
  }

  e_coef <- e_se <- array(0, c(reps, ntaus)) # estimated with QR
  s_coef <- s_se <- array(0, c(reps, ntaus)) # smoothed est.
  c_coef <- c_se <- array(0, c(reps, ntaus)) # conquer
  m_coef <- array(0, c(reps, nboot, ntaus))
  m_se <- array(0, c(reps, ntaus))

  for (j in 1:reps) {
    dat <- dgp(n, T, nml = nml, het = het, dep = dep, dyn = dyn)
    y <- dat$y
    x <- dat$x
    w <- rexp(n * nboot, 1)
    U <- matrix(rep(w, each = T), n * T, nboot) # assumes T_i = T for all i

    for (k in 1:ntaus) {
      rqest <- rq(y ~ x + d - 1, tau = taus[k])
      sm <- srq(y, x, d, rqest, n, T, taus[k])
      cq <- conq(y, x, d, n, T, taus[k])
      e_coef[j, k] <- rqest$coef[1]
      e_se[j, k] <- est_se(rqest$coef[1], rqest$resid, x, n, T, tau = taus[k])
      s_coef[j, k] <- sm$co[1]
      s_se[j, k] <- est_se(sm$coef[1], sm$resid, x, n, T, tau = taus[k])
      c_coef[j, k] <- cq$co[1]
      c_se[j, k] <- est_se(cq$coef, cq$resid, x, n, T, tau = taus[k])
      m_coef[j, , k] <- boot.rq.wxy(cbind(x, d), y, U, tau = taus[k])[, 1]
      m_se[j, k] <- sd(m_coef[j, , k])
    }
  }

  sim_ans <- list(e_coef = e_coef, s_coef = s_coef, c_coef = c_coef,
                  e_se = e_se, s_se = s_se, c_se = c_se, m_se = m_se,
                  taus = taus, true_beta = true_beta, n = n, T = T,
                  alpha = alpha, nml = nml, het = het, dep = dep, dyn = dyn)
  sim_ans
}

.lec.CurrentStream(stream_names[pnum])
sims <- sim_run(100, 100, dnum, simreps, 2, taus, al)
.lec.CurrentStreamEnd()

save(sims, file = paste0("./parts/design", dnum, "part", pnum, ".rda"))

