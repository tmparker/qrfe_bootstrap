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

# Run a simulation experiment.  des_num is the design number.
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

  # Several estimators here
  e_coef <- e_se <- e_cover <- array(0, c(reps, ntaus)) # estimated with QR
  s_coef <- s_se <- s_cover <- array(0, c(reps, ntaus)) # smoothed est.
  c_coef <- c_se <- c_cover <- array(0, c(reps, ntaus)) # conquer
  ebc_coef <- ebc_se <- ebc_cover <- array(0, c(reps, ntaus)) # bias-corrected rq
  sbc_coef <- sbc_se <- sbc_cover <- array(0, c(reps, ntaus)) # bc smoothed rq
  j_coef <- j_cover <- array(0, c(reps, ntaus))
  m_coef <- m_t <- array(0, c(reps, nboot, ntaus))
  m_mean <- m_se <- array(0, c(reps, ntaus))
  m_pCI <- m_seCI <- m_tCI <- array(0, c(reps, 2, ntaus))
  m_cover_p <- m_cover_se <- m_cover_t <- array(0, c(reps, ntaus))

  for (j in 1:reps) {
    dat <- dgp(n, T, nml = nml, het = het, dep = dep, dyn = dyn)
    y <- dat$y
    x <- dat$x
    w <- rexp(n * nboot, 1)
    U <- matrix(rep(w, each = T), n * T, nboot) # assumes T_i = T for all i

    for (k in 1:ntaus) {
      rqest <- rq.fit(cbind(x, d), y, tau = taus[k])
      e_coef[j, k] <- rqest$coef[1]
      e_se[j, k] <- est_se(e_coef[j, k], rqest$resid, x, n, T, tau = taus[k])
      e_cover[j, k] <- 1 * (true_beta[k] > e_coef[j, k] - z_al2 * e_se[j, k] &
                            true_beta[k] < e_coef[j, k] + z_al2 * e_se[j, k])
      # Smoothed estimators
      # Galvao & Kato 2016 smoothed est.
      sm <- srq(y, x, d, rqest, n, T, taus[k])
      s_coef[j, k] <- sm$coef[1]
      s_se[j, k] <- est_se(s_coef[j, k], sm$resid, x, n, T, tau = taus[k])
      s_cover[j, k] <- 1 * (true_beta[k] > s_coef[j, k] - z_al2 * s_se[j, k] &
                            true_beta[k] < s_coef[j, k] + z_al2 * s_se[j, k])

      # Bias-corrected estimates (only smoothed is trustworthy)
      bcest <- srq_biascorrect(y, x, d, n, T, taus[k])
      ebc_coef[j, k] <- bcest[1]
      ebc_cover[j, k] <- 1 * (true_beta[k] > ebc_coef[j, k] - z_al2 * e_se[j, k] &
                            true_beta[k] < ebc_coef[j, k] + z_al2 * e_se[j, k])
      sbc_coef[j, k] <- bcest[2]
      sbc_cover[j, k] <- 1 * (true_beta[k] > sbc_coef[j, k] - z_al2 * s_se[j, k] &
                            true_beta[k] < sbc_coef[j, k] + z_al2 * s_se[j, k])
      # Jackknife bias-corrected?  Use smoothed reg. standard error est.
      j_coef[j, k] <- half_panel_jackknife(y, x, d, n, T, taus[k])
      j_cover[j, k] <- 1 * (true_beta[k] > j_coef[j, k] - z_al2 * s_se[j, k] &
                            true_beta[k] < j_coef[j, k] + z_al2 * s_se[j, k])

      # Multiplier bootstrap
      m_coef[j, , k] <- boot.rq.wxy(cbind(x, d), y, U, tau = taus[k])[, 1]
      m_mean[j, k] <- mean(m_coef[j, , k])
      m_se[j, k] <- sd(m_coef[j, , k])
      # bootstrap t statistic
      m_t[j, , k] <- (m_coef[j, , k] - e_coef[j, k]) / m_se[j, k]
      # percentile CI
      m_pCI[j, , k] <- c(quantile(m_coef[j, , k], al2),
                          quantile(m_coef[j, , k], 1 - al2))
      m_cover_p[j, k] <- 1 * (true_beta[k] > m_pCI[j, 1, k] &
                              true_beta[k] < m_pCI[j, 2, k])
      # bootstrap standard error estimation, normal quantile
      m_seCI[j, , k] <- c(e_coef[j, k] - z_al2 * m_se[j, k],
                          e_coef[j, k] + z_al2 * m_se[j, k])
      m_cover_se[j, k] <- 1 * (true_beta[k] > m_seCI[j, 1, k] &
                                true_beta[k] < m_seCI[j, 2, k])
      # bootstrap reference distribution, estimated s.e.
      m_tCI[j, , k] <- c(e_coef[j, k] - quantile(m_t[j, , k], 1 - al2) * e_se[j, k],
                          e_coef[j, k] - quantile(m_t[j, , k], al2) * e_se[j, k])
      m_cover_t[j, k] <- 1 * (true_beta[k] > m_tCI[j, 1, k] &
                              true_beta[k] < m_tCI[j, 2, k])

    }
  }
  sim_ans <- list(e_coef = e_coef, e_se = e_se, e_cover = e_cover,
                  s_coef = s_coef, s_se = s_se, s_cover = s_cover,
                  ebc_coef = ebc_coef, ebc_cover = ebc_cover,
                  sbc_coef = sbc_coef, sbc_cover = sbc_cover,
                  j_coef = j_coef, j_cover = j_cover,
                  m_coef = m_coef, m_t = m_t, m_mean = m_mean, m_se = m_se,
                  m_pCI = m_pCI, m_seCI = m_seCI, m_tCI = m_tCI,
                  m_cover_p = m_cover_p, m_cover_se = m_cover_se,
                  m_cover_t = m_cover_t,
                  taus = taus, true_beta = true_beta,
                  n = n, T = T, alpha = alpha, nml = nml,
                  het = het, dep = dep, dyn = dyn)
  sim_ans
}

.lec.CurrentStream(stream_names[pnum])
sims <- sim_run(100, 100, dnum, simreps, nboot, taus, al)
.lec.CurrentStreamEnd()

save(sims, file = paste0("./parts/design", dnum, "part", pnum, ".rda"))

