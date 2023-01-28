# This script generates data and uses the bootstrap to conduct inference about
# the parameters.  Several types of measurements are made to facilitate
# comparisons later.  The actual comparisons are processed in another script.

rm(list=ls())
library(quantreg)
library(rlecuyer)

source("../../common.R")

# 100 instances, each instance runs 10 repetitions.
# Simulation design 1 means arg should run 0-99
arg <- as.numeric(commandArgs(trailingOnly = TRUE))
nstream <- 100 # number of parallel random number streams
dnum <- ceiling((arg + 1) / nstream) # design number
pnum <- arg %% nstream + 1
stream_names <- paste0("st", 1:nstream)
.lec.SetPackageSeed(6:1)
.lec.CreateStream(stream_names) # Seed info stored in list .lec.Random.seed.table

# Run a simulation experiment.  des_num is the design number.
sim_run <- function(n, T, des_num, reps, nboot, taus = 1:3 / 4, alpha = 0.1) {
  ntaus <- length(taus)
  al2 <- alpha / 2
  z_al2 <- qnorm(1 - al2)
  d <- kronecker(diag(n), rep(1, T))
  est <- vector(mode = "list", length = ntaus)

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

  e_coef <- e_se <- array(0, c(reps, ntaus))
  e_cover <- array(0, c(reps, ntaus))
  m_coef <- m_t <- array(0, c(reps, nboot, ntaus))
  m_mean <- m_se <- array(0, c(reps, ntaus))
  m_pCI1 <- m_pCI2 <- array(0, c(reps, 2, ntaus))
  m_cover_p1 <- m_cover_p2 <- array(0, c(reps, ntaus))

  for (j in 1:reps) {
    dat <- dgp(n, T, nml = nml, het = het, dep = dep, dyn = dyn)
    y <- dat$y
    x <- dat$x
    # bootstrap weights generated up front for all the bootstrap repetitions
    w <- rexp(n * nboot, 1)
    U <- matrix(rep(w, each = T), n * T, nboot) # assumes T_i = T for all i

    for (k in 1:ntaus) {
      rqest <- rq(y ~ x + d - 1, tau = taus[k])
      # e_ for estimate (asymptotic)
      e_coef[j, k] <- rqest$coef[1]
      e_se[j, k] <- est_se(rqest$coef[1], rqest$resid, x, n, T, tau = taus[k])
      e_cover[j, k] <- 1 * (true_beta[k] > e_coef[j, k] - z_al2 * e_se[j, k] &
                            true_beta[k] < e_coef[j, k] + z_al2 * e_se[j, k])
      # m_ for multiplier (bootstrap)
      m_coef[j, , k] <- boot.rq.wxy(cbind(x, d), y, U, tau = taus[k])[, 1]
      m_mean[j, k] <- mean(m_coef[j, , k])
      m_se[j, k] <- sd(m_coef[j, , k])
      m_t[j, , k] <- (m_coef[j, , k] - e_coef[j, k]) / m_se[j, k]
      m_pCI1[j, , k] <- c(2 * e_coef[j, k] - quantile(m_coef[j, , k], 1 - al2),
                          2 * e_coef[j, k] - quantile(m_coef[j, , k], al2))
      m_pCI2[j, , k] <- c(quantile(m_coef[j, , k], al2),
                          quantile(m_coef[j, , k], 1 - al2))
      m_cover_p1[j, k] <- 1 * (true_beta[k] > m_pCI1[j, 1, k] &
                              true_beta[k] < m_pCI1[j, 2, k])
      m_cover_p2[j, k] <- 1 * (true_beta[k] > m_pCI2[j, 1, k] &
                              true_beta[k] < m_pCI2[j, 2, k])
    }
  }

  sim_ans <- list(e_coef = e_coef, e_se = e_se, e_cover = e_cover,
                  m_coef = m_coef, m_t = m_t, m_mean = m_mean, m_se = m_se,
                  m_pCI1 = m_pCI1, m_pCI2 = m_pCI2, m_cover_p1 = m_cover_p1,
                  m_cover_p2 = m_cover_p2, taus = taus, true_beta = true_beta,
                  n = n, T = T, alpha = alpha, nml = nml, het = het, dep = dep,
                  dyn = dyn)
  sim_ans
}

.lec.CurrentStream(stream_names[pnum])
sims <- sim_run(100, 100, dnum, simreps, nboot, taus, al)
.lec.CurrentStreamEnd()

save(sims, file = paste0("./parts/design", dnum, "part", pnum, ".rda"))

