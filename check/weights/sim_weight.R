# This is focused on bootstrap weight distributions.  The distributions are
# exponential, chi^2, uniform and folded normal.  The last three of these
# choices require more careful tuning than the exponential, but can be used.
# This script generates data, estimates models and conducts inference over
# simulation repetitions.  Comparisons are processed in another script.

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

cw_df <- 5 # df of the chi-square weights - arbitrary

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

  e_coef <- e_se <- array(0, c(reps, ntaus))
  e_cover <- array(0, c(reps, ntaus))

  m_coef <- m_t <- array(0, c(reps, nboot, ntaus)) # exponential weights
  m_mean <- m_se <- array(0, c(reps, ntaus))
  m_pCI <- m_pCI2 <- m_seCI <- m_tCI <- array(0, c(reps, 2, ntaus))
  m_cover_p <- m_cover_p2 <- m_cover_se <- m_cover_t <- array(0, c(reps, ntaus))

  c_coef <- c_t <- array(0, c(reps, nboot, ntaus)) # chi square / df
  c_mean <- c_se <- array(0, c(reps, ntaus))
  c_pCI <- c_pCI2 <- c_seCI <- c_tCI <- array(0, c(reps, 2, ntaus))
  c_cover_p <- c_cover_p2 <- c_cover_se <- c_cover_t <- array(0, c(reps, ntaus))
  sc <- sqrt(2 / cw_df) # sd of chi square weights (w = X / df)

  u_coef <- u_t <- array(0, c(reps, nboot, ntaus)) # uniform
  u_mean <- u_se <- array(0, c(reps, ntaus))
  u_pCI <- u_pCI2 <- u_seCI <- u_tCI <- array(0, c(reps, 2, ntaus))
  u_cover_p <- u_cover_p2 <- u_cover_se <- u_cover_t <- array(0, c(reps, ntaus))
  su <- sqrt(1 / 3) # sd of U[0, 2] weights

  n_coef <- n_t <- array(0, c(reps, nboot, ntaus)) # folded normal
  n_mean <- n_se <- array(0, c(reps, ntaus))
  n_pCI <- n_pCI2 <- n_seCI <- n_tCI <- array(0, c(reps, 2, ntaus))
  n_cover_p <- n_cover_p2 <- n_cover_se <- n_cover_t <- array(0, c(reps, ntaus))
  sn <- sqrt(pi / 2 - 1) # sd of folded normal (s.t. E[W] = 1)

  for (j in 1:reps) {
    dat <- dgp(n, T, nml = nml, het = het, dep = dep, dyn = dyn)
    y <- dat$y
    x <- dat$x
    w0 <- rexp(n * nboot, 1)
    w1 <- rchisq(n * nboot, df = cw_df) / cw_df # so mean = 1
    w2 <- runif(n * nboot, min = 0, max = 2)
    w3 <- abs(rnorm(n * nboot, sd = sqrt(pi / 2))) # also so mean = 1
    E <- matrix(rep(w0, each = T), n * T, nboot)
    Ch <- matrix(rep(w1, each = T), n * T, nboot)
    U <- matrix(rep(w2, each = T), n * T, nboot)
    FN <- matrix(rep(w3, each = T), n * T, nboot)

    for (k in 1:ntaus) {
      rqest <- rq(y ~ x + d - 1, tau = taus[k])
      e_coef[j, k] <- rqest$coef[1]
      e_se[j, k] <- est_se(rqest$coef[1], rqest$resid, x, n, T, tau = taus[k])
      e_cover[j, k] <- 1 * (true_beta[k] > e_coef[j, k] - z_al2 * e_se[j, k] &
                            true_beta[k] < e_coef[j, k] + z_al2 * e_se[j, k])

      # Exponential
      m_coef[j, , k] <- boot.rq.wxy(cbind(x, d), y, E, tau = taus[k])[, 1]
      m_mean[j, k] <- mean(m_coef[j, , k])
      m_se[j, k] <- sd(m_coef[j, , k])
      # bootstrap t statistic
      m_t[j, , k] <- (m_coef[j, , k] - e_coef[j, k]) / m_se[j, k]
      # percentile CI
      m_pCI[j, , k] <- c(quantile(m_coef[j, , k], al2),
                          quantile(m_coef[j, , k], 1 - al2))
      m_cover_p[j, k] <- 1 * (true_beta[k] > m_pCI[j, 1, k] &
                              true_beta[k] < m_pCI[j, 2, k])
      # Test the other more "proper" CI one more time
      m_pCI2[j, , k] <- c(2 * e_coef[j, k] - quantile(m_coef[j, , k], 1 - al2),
                          2 * e_coef[j, k] - quantile(m_coef[j, , k], al2))
      m_cover_p2[j, k] <- 1 * (true_beta[k] > m_pCI2[j, 1, k] &
                              true_beta[k] < m_pCI2[j, 2, k])
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

      # w = chi square / its df
      c_coef[j, , k] <- boot.rq.wxy(cbind(x, d), y, Ch, tau = taus[k])[, 1]
      c_mean[j, k] <- mean(c_coef[j, , k])
      c_se[j, k] <- sd(c_coef[j, , k])
      # bootstrap t statistic
      c_t[j, , k] <- (c_coef[j, , k] - e_coef[j, k]) / c_se[j, k]
      # percentile CI (all other dbns different from exponential)
      c_pCI[j, , k] <- c((1 - 1 / sc) * e_coef[j, k]
                          + quantile(c_coef[j, , k], al2) / sc,
                          (1 - 1 / sc) * e_coef[j, k]
                          + quantile(c_coef[j, , k], 1 - al2) / sc)
      c_cover_p[j, k] <- 1 * (true_beta[k] > c_pCI[j, 1, k] &
                              true_beta[k] < c_pCI[j, 2, k])
      c_pCI2[j, , k] <- c((1 + 1 / sc) * e_coef[j, k]
                          - quantile(c_coef[j, , k], 1 - al2) / sc,
                          (1 + 1 / sc) * e_coef[j, k]
                          - quantile(c_coef[j, , k], al2) / sc)
      c_cover_p2[j, k] <- 1 * (true_beta[k] > c_pCI2[j, 1, k] &
                              true_beta[k] < c_pCI2[j, 2, k])
      # bootstrap standard error estimation, normal quantile
      c_seCI[j, , k] <- c(e_coef[j, k] - z_al2 * c_se[j, k] / sc,
                          e_coef[j, k] + z_al2 * c_se[j, k] / sc)
      c_cover_se[j, k] <- 1 * (true_beta[k] > c_seCI[j, 1, k] &
                                true_beta[k] < c_seCI[j, 2, k])
      # bootstrap reference distribution, estimated s.e.
      c_tCI[j, , k] <- c(e_coef[j, k] - quantile(c_t[j, , k], 1 - al2) * e_se[j, k],
                          e_coef[j, k] - quantile(c_t[j, , k], al2) * e_se[j, k])
      c_cover_t[j, k] <- 1 * (true_beta[k] > c_tCI[j, 1, k] &
                              true_beta[k] < c_tCI[j, 2, k])

      # U[0, 2]
      u_coef[j, , k] <- boot.rq.wxy(cbind(x, d), y, U, tau = taus[k])[, 1]
      u_mean[j, k] <- mean(u_coef[j, , k])
      u_se[j, k] <- sd(u_coef[j, , k])
      # bootstrap t statistic
      u_t[j, , k] <- (u_coef[j, , k] - e_coef[j, k]) / u_se[j, k]
      # percentile CI
      u_pCI[j, , k] <- c((1 - 1 / su) * e_coef[j, k]
                          + quantile(u_coef[j, , k], al2) / su,
                          (1 - 1 / su) * e_coef[j, k]
                          + quantile(u_coef[j, , k], 1 - al2) / su)
      u_cover_p[j, k] <- 1 * (true_beta[k] > u_pCI[j, 1, k] &
                              true_beta[k] < u_pCI[j, 2, k])
      u_pCI2[j, , k] <- c((1 + 1 / su) * e_coef[j, k]
                          - quantile(u_coef[j, , k], 1 - al2) / su,
                          (1 + 1 / su) * e_coef[j, k]
                          - quantile(u_coef[j, , k], al2) / su)
      u_cover_p2[j, k] <- 1 * (true_beta[k] > u_pCI2[j, 1, k] &
                              true_beta[k] < u_pCI2[j, 2, k])
      # bootstrap standard error estimation, normal quantile
      u_seCI[j, , k] <- c(e_coef[j, k] - z_al2 * u_se[j, k] / su,
                          e_coef[j, k] + z_al2 * u_se[j, k] / su)
      u_cover_se[j, k] <- 1 * (true_beta[k] > u_seCI[j, 1, k] &
                                true_beta[k] < u_seCI[j, 2, k])
      # bootstrap reference distribution, estimated s.e.
      u_tCI[j, , k] <- c(e_coef[j, k] - quantile(u_t[j, , k], 1 - al2) * e_se[j, k],
                          e_coef[j, k] - quantile(u_t[j, , k], al2) * e_se[j, k])
      u_cover_t[j, k] <- 1 * (true_beta[k] > u_tCI[j, 1, k] &
                              true_beta[k] < u_tCI[j, 2, k])

      # Folded normal
      n_coef[j, , k] <- boot.rq.wxy(cbind(x, d), y, FN, tau = taus[k])[, 1]
      n_mean[j, k] <- mean(n_coef[j, , k])
      n_se[j, k] <- sd(n_coef[j, , k])
      # bootstrap t statistic
      n_t[j, , k] <- (n_coef[j, , k] - e_coef[j, k]) / n_se[j, k]
      # percentile CI
      n_pCI[j, , k] <- c((1 - 1 / sn) * e_coef[j, k]
                          + quantile(n_coef[j, , k], al2) / sn,
                          (1 - 1 / sn) * e_coef[j, k]
                          + quantile(n_coef[j, , k], 1 - al2) / sn)
      n_cover_p[j, k] <- 1 * (true_beta[k] > n_pCI[j, 1, k] &
                              true_beta[k] < n_pCI[j, 2, k])
      n_pCI2[j, , k] <- c((1 + 1 / sn) * e_coef[j, k]
                          - quantile(n_coef[j, , k], 1 - al2) / sn,
                          (1 + 1 / sn) * e_coef[j, k]
                          - quantile(n_coef[j, , k], al2) / sn)
      n_cover_p2[j, k] <- 1 * (true_beta[k] > n_pCI2[j, 1, k] &
                              true_beta[k] < n_pCI2[j, 2, k])
      # bootstrap standard error estimation, normal quantile
      n_seCI[j, , k] <- c(e_coef[j, k] - z_al2 * n_se[j, k] / sn,
                          e_coef[j, k] + z_al2 * n_se[j, k] / sn)
      n_cover_se[j, k] <- 1 * (true_beta[k] > n_seCI[j, 1, k] &
                                true_beta[k] < n_seCI[j, 2, k])
      # bootstrap reference distribution, estimated s.e.
      n_tCI[j, , k] <- c(e_coef[j, k] - quantile(n_t[j, , k], 1 - al2) * e_se[j, k],
                          e_coef[j, k] - quantile(n_t[j, , k], al2) * e_se[j, k])
      n_cover_t[j, k] <- 1 * (true_beta[k] > n_tCI[j, 1, k] &
                              true_beta[k] < n_tCI[j, 2, k])

    }
  }
  sim_ans <- list(e_coef = e_coef, e_se = e_se, e_cover = e_cover,
                  m_coef = m_coef, m_t = m_t, m_mean = m_mean, m_se = m_se,
                  m_pCI = m_pCI, m_pCI2 = m_pCI2, m_seCI = m_seCI, m_tCI = m_tCI,
                  m_cover_p = m_cover_p, m_cover_p2 = m_cover_p2,
                  m_cover_se = m_cover_se, m_cover_t = m_cover_t,

                  c_coef = c_coef, c_t = c_t, c_mean = c_mean, c_se = c_se,
                  c_pCI = c_pCI, c_seCI = c_seCI, c_tCI = c_tCI,
                  c_cover_p = c_cover_p, c_cover_p2 = c_cover_p2,
                  c_cover_se = c_cover_se, c_cover_t = c_cover_t,

                  u_coef = u_coef, u_t = u_t, u_mean = u_mean, u_se = u_se,
                  u_pCI = u_pCI, u_seCI = u_seCI, u_tCI = u_tCI,
                  u_cover_p = u_cover_p, u_cover_p2 = u_cover_p2,
                  u_cover_se = u_cover_se, u_cover_t = u_cover_t,

                  n_coef = n_coef, n_t = n_t, n_mean = n_mean, n_se = n_se,
                  n_pCI = n_pCI, n_seCI = n_seCI, n_tCI = n_tCI,
                  n_cover_p = n_cover_p, n_cover_p2 = n_cover_p2,
                  n_cover_se = n_cover_se, n_cover_t = n_cover_t,

                  taus = taus, true_beta = true_beta,
                  n = n, T = T, alpha = alpha, nml = nml,
                  het = het, dep = dep, dyn = dyn)
                  print(paste("n = ", n, "T = ", T))
  sim_ans
}

.lec.CurrentStream(stream_names[pnum])
sims <- sim_run(100, 100, dnum, simreps, nboot, taus, al)
.lec.CurrentStreamEnd()

save(sims, file = paste0("./parts/design", dnum, "part", pnum, ".rda"))

