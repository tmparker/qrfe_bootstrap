# This script generates data and uses the bootstrap to conduct inference about
# the parameters.  Several types of measurements are made to facilitate
# comparisons later.  The actual comparisons are processed in another script.
# This one just provides the raw data.

rm(list=ls())
library(quantreg)
library(rlecuyer)

# dnum and pnum specify the design and part numbers.
# Designs are numbered as described below.
# arg should go from 100 to 1099 - 10 designs, 100 subparts per design
arg <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
dnum <- floor(arg / 100)
pnum <- arg %% 100 + 1

nstream <- 100 # number of parallel random number streams
stream_names <- paste0("st", 1:nstream)
.lec.SetPackageSeed(6:1)
.lec.CreateStream(stream_names) # Seed info stored in list .lec.Random.seed.table
# 100 different per-part RNG seeds set below around the simulation run.

# These parameters are constant across all the simulations
taus <- 1:3 / 4
alpha <- 0.1
simreps <- 10

nboot <- 999
nvec <- c(25, 50, 100)
tvec <- c(20, 50, 100)

# List of designs
# Designs are numbered 1-10, baed on alphabetical order of names.  For example,
# chisq_dep_het is the first one alphabetically, so it's design 1.  The result
# of this function is used to set the arguments for the subsequent functions.
# Unfortunately the naming convention is different for the file labels and the
# order of the arguments here: heteroskedasticity comes before dependence in the
# arguments, but in the file names, dependence comes first (only matters for
# organizing tables).
design <- function(k) {
  nml <- FALSE
  het <- FALSE
  dep <- FALSE
  dyn <- FALSE
  if (k %in% 6:10) nml <- TRUE
  if (k %in% c(1, 4, 6, 9)) het <- TRUE
  if (k %in% c(1, 2, 6, 7)) dep <- TRUE
  if (k %in% c(3, 8)) dyn <- TRUE
  return(list(nml = nml, het = het, dep = dep, dyn = dyn))
}

# Generate a sample.
# By default uses iid normals.  Change to chi^2, heteroskedastic, dependent
# errors by changing the logical arguments.
dgp <- function(n, t, nml = TRUE, het = FALSE, dep = FALSE, dyn = FALSE,
                be = 1, ga = 0.2, the = 0.5, rho = 0.4, dyn_burnin = 50) {
  if (dyn) {
    al <- runif(n)
    u <- matrix(rnorm(n * (t + dyn_burnin)), t + dyn_burnin, n)
    yaux <- matrix(0, t + dyn_burnin, n)
    yaux[1, ] <- al + u[1, ]
    for(i in 2:(t + dyn_burnin)) {
      yaux[i, ] <- al + rho * yaux[i-1, ] + u[i, ]
    }
    yaux.y <- yaux[(dyn_burnin + 1):(t + dyn_burnin), ]
    yaux.yl <- yaux[dyn_burnin:(t + dyn_burnin-1), ]
    y <- c(yaux.y)
    x <- c(yaux.yl)
    return(list(y = y, x = x, true_al = al))
  } else {
    if (!het) { ga <- 0 }
    al <- rep(runif(n), each = t)
    x <- rchisq(n * t, df = 3) + 0.3 * al
    if (dep) { # MA(1) error terms if dependent:
      if (nml) {
        e <- rnorm(n * (t + 1))
      } else {
        e <- rchisq(n * (t + 1), df = 4)
      }
      i1 <- rep((1:n - 1) * (t + 1), each = t) + 2:(t + 1)
      i0 <- rep((1:n - 1) * (t + 1), each = t) + 1:t
      u <- e[i1] + the * e[i0]
    } else { # for independent errors
      if (nml) {
        u <- rnorm(n * t)
      } else {
        u <- rchisq(n * t, df = 4)
      }
    }
    y <- al + be * x + (1 + ga * x) * u
    return(list(y = y, x = x, true_al = al))
  }
}

# The true parameters for the dynamic specification we use are the same as for
# the locations-shift specification, so there's no option for a dynamic model.
# Just use beta_true(taus, be = 0.4), since that is the default value for the
# dynamic DGP.
beta_true <- function(taus, nml = TRUE, het = FALSE, dep = FALSE,
                      be = 1, ga = 0.2, the = 0.5) {
  if (nml) { Q <- qnorm(taus) } else { Q <- qchisq(taus, df = 4) }
  if (!het) { ga <- 0 }
  if (dep) {
    if (nml) {
      Q <- qnorm(taus, sd = sqrt(1 + the^2)) # easy if normal
    } else {
      e <- rchisq(1000001, df = 4) # simulate if chi^2 MA(1)
      u <- e[-1] + the * e[-length(e)]
      Q <- quantile(u, probs = taus)
    }
  }
  btrue <- be + ga * Q
  btrue
}

# For checking the accuracy of standard error estimates.  Just get a
# monte carlo approximation of the true standard error of betahat rather than
# using the asymptotic formula.  Should be done conditional on x and intercepts,
# but al has to be the estimated alpha so that there is a way to generate ys.
# For the dynamic model, can't make it conditional on lagged responses
# observations, but it can still be conditional on (estimated) alpha.
se_true <- function(n, t, taus, x, al, nml = TRUE, het = FALSE, dep = FALSE,
                    dyn = FALSE, be = 1, ga = 0.2, the = 0.5, rho = 0.4,
                    dyn_burnin = 50, sim_reps = 1000) {
  D <- kronecker(diag(n), rep(1, t))
  ests <- matrix(0, sim_reps, length(taus))
  if (dyn) {
    for (r in 1:sim_reps) {
      u <- matrix(rnorm(n * (t + dyn_burnin)), t + dyn_burnin, n)
      yaux <- matrix(0, t + dyn_burnin, n)
      # Assumed that alpha is like runif(n)
      yaux[1, ] <- al + u[1, ]
      for(i in 2:(t + dyn_burnin)) {
        yaux[i, ] <- al + rho * yaux[i-1, ] + u[i, ]
      }
      yaux.y <- yaux[(dyn_burnin + 1):(t + dyn_burnin), ]
      yaux.yl <- yaux[dyn_burnin:(t + dyn_burnin-1), ]
      y <- c(yaux.y)
      x <- c(yaux.yl)
      ests[r, ] <- rq(y ~ x + D - 1, tau = taus)$coef["x", drop = FALSE]
    }
  } else { # all the not-dynamic DGPs here
    for (r in 1:sim_reps) {
      if (!het) { ga <- 0 }
      if (dep) { # MA(1) error terms if dependent:
        if (nml) {
          e <- rnorm(n * (t + 1))
        } else {
          e <- rchisq(n * (t + 1), df = 4)
        }
        i1 <- rep((1:n - 1) * (t + 1), each = t) + 2:(t + 1)
        i0 <- rep((1:n - 1) * (t + 1), each = t) + 1:t
        u <- e[i1] + the * e[i0]
      } else { # for independent errors
        if (nml) {
          u <- rnorm(n * t)
        } else {
          u <- rchisq(n * t, df = 4)
        }
      }
      # Assumed that al is like rep(runif(n), each = t).
      y <- al + be * x + (1 + ga * x) * u
      ests[r, ] <- rq(y ~ x + D - 1, tau = taus)$coef["x", drop = FALSE]
    }
  }
  true_se <- apply(ests, 2, sd) # one estimate for each quantile
  true_se
}

# Run a simulation experiment.  des_num is the design number.
sim_run <- function(n, t, des_num, reps, nboot, taus = 1:3 / 4, alpha = 0.1) {
  ntaus <- length(taus)
  al2 <- alpha / 2
  est <- vector(mode = "list", length = ntaus)
  z_al2 <- qnorm(1 - al2)
  D <- kronecker(diag(n), rep(1, t))

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

  # We've got all kinds of results!  They will be processed later.
  e_coef <- e_sd <- array(0, c(reps, ntaus))
  e_cover <- array(0, c(reps, ntaus))
  m_coef <- m_t <- array(0, c(reps, nboot, ntaus))
  m_mean <- m_sd <- array(0, c(reps, ntaus))
  m_pCI <- m_seCI <- m_tCI <- array(0, c(reps, 2, ntaus))
  m_cover_p <- m_cover_se <- m_cover_t <- array(0, c(reps, ntaus))
  g_coef <- g_t <- array(0, c(reps, nboot, ntaus))
  g_mean <- g_sd <- array(0, c(reps, ntaus))
  g_pCI <- g_seCI <- g_tCI <- array(0, c(reps, 2, ntaus))
  g_cover_p <- g_cover_se <- g_cover_t <- array(0, c(reps, ntaus))

  for (j in 1:reps) {
    dat <- dgp(n, t, nml = nml, het = het, dep = dep, dyn = dyn)
    y <- dat$y
    x <- dat$x
    # "True standard errors" found using Monte Carlo
    true_se <- se_true(n, t, taus, x, dat$true_al, nml = nml, het = het,
                        dep = dep, dyn = dyn)
    # bootstrap weights generated up front for all the bootstrap repetitions
    w <- rexp(n * nboot, 1)
    U <- matrix(rep(w, each = t), n * t, nboot) # assumes T_i = T for all i

    for (k in 1:ntaus) {
      # e_xxx stands for "estimate", as in the estimated coefficient
      est[[k]] <- summary(rq(y ~ x + D - 1, tau = taus[k]), se = "ker")
      e_coef[j, k] <- est[[k]]$coef[1, 1]
      e_sd[j, k] <- est[[k]]$coef[1, 2]
      e_cover[j, k] <- 1 * (true_beta[k] > e_coef[j, k] - z_al2 * e_sd[j, k] &
                            true_beta[k] < e_coef[j, k] + z_al2 * e_sd[j, k])

      # m_xxx stands for multiplier bootstrap applied to the objective function
      m_coef[j, , k] <- boot.rq.wxy(cbind(x, D), y, U, tau = taus[k])[, 1]
      m_mean[j, k] <- mean(m_coef[j, , k])
      m_sd[j, k] <- sd(m_coef[j, , k])
      # using a bootstrap estimate of the standard error of betastar
      m_t[j, , k] <- (m_coef[j, , k] - e_coef[j, k]) / m_sd[j, k]
      # Bootstrap percentile intervals_
      m_pCI[j, , k] <- c(quantile(m_coef[j, , k], al2),
                          quantile(m_coef[j, , k], 1 - al2))
      m_cover_p[j, k] <- 1 * (true_beta[k] > m_pCI[j, 1, k] &
                              true_beta[k] < m_pCI[j, 2, k])
      # Use bootstrap for standard error estimation, normal quantile for endpoints
      m_seCI[j, , k] <- c(e_coef[j, k] - z_al2 * m_sd[j, k],
                          e_coef[j, k] + z_al2 * m_sd[j, k])
      m_cover_se[j, k] <- 1 * (true_beta[k] > m_seCI[j, 1, k] &
                                true_beta[k] < m_seCI[j, 2, k])
      # Use bootstrap for reference distribution, estimated s.e.
      m_tCI[j, , k] <- c(e_coef[j, k] - quantile(m_t[j, , k], 1 - al2) * e_sd[j, k],
                          e_coef[j, k] - quantile(m_t[j, , k], al2) * e_sd[j, k])
      m_cover_t[j, k] <- 1 * (true_beta[k] > m_tCI[j, 1, k] &
                              true_beta[k] < m_tCI[j, 2, k])

      # g_xxx stands for gradient bootstrap
      g_coef[j, , k] <- boot.rq(cbind(x, D), y, tau = taus[k], R = nboot,
                                cluster = rep(1:n, rep(t, n)))$B[,1]
      g_mean[j, k] <- mean(g_coef[j, , k])
      g_sd[j, k] <- sd(g_coef[j, , k])
      # using a bootstrap estimate of the standard error of betastar!
      g_t[j, , k] <- (g_coef[j, , k] - e_coef[j, k]) / g_sd[j, k]
      # Bootstrap percentile intervals_
      g_pCI[j, , k] <- c(quantile(g_coef[j, , k], al2),
                          quantile(g_coef[j, , k], 1 - al2))
      g_cover_p[j, k] <- 1 * (true_beta[k] > g_pCI[j, 1, k] &
                              true_beta[k] < g_pCI[j, 2, k])
      # Use bootstrap for standard error estimation, normal quantile for endpoints
      g_seCI[j, , k] <- c(e_coef[j, k] - z_al2 * g_sd[j, k],
                          e_coef[j, k] + z_al2 * g_sd[j, k])
      g_cover_se[j, k] <- 1 * (true_beta[k] > g_seCI[j, 1, k] &
                                true_beta[k] < g_seCI[j, 2, k])
      # Use bootstrap for reference distribution, estimated s.e.
      g_tCI[j, , k] <- c(e_coef[j, k] - quantile(g_t[j, , k], 1 - al2) * e_sd[j, k],
                          e_coef[j, k] - quantile(g_t[j, , k], al2) * e_sd[j, k])
      g_cover_t[j, k] <- 1 * (true_beta[k] > g_tCI[j, 1, k] &
                              true_beta[k] < g_tCI[j, 2, k])
    }
  }

  sim_ans <- list(e_coef = e_coef, e_sd = e_sd, e_cover = e_cover,
                  m_coef = m_coef, m_t = m_t, m_mean = m_mean, m_sd = m_sd,
                  m_pCI = m_pCI, m_seCI = m_seCI, m_tCI = m_tCI,
                  m_cover_p = m_cover_p, m_cover_se = m_cover_se,
                  m_cover_t = m_cover_t, g_coef = g_coef, g_t = g_t,
                  g_mean = g_mean, g_sd = g_sd, g_pCI = g_pCI, g_seCI = g_seCI,
                  g_tCI = g_tCI, g_cover_p = g_cover_p, g_cover_se = g_cover_se,
                  g_cover_t = g_cover_t, taus = taus, true_beta = true_beta,
                  true_se = true_se, n = n, t = t, alpha = alpha, nml = nml,
                  het = het, dep = dep, dyn = dyn)
  sim_ans
}

.lec.CurrentStream(stream_names[pnum])
sims <- mapply(FUN = sim_run, rep(nvec, each = length(tvec)),
                rep(tvec, length(nvec)), MoreArgs = list(des_num = dnum,
                reps = simreps, nboot = nboot, taus = taus, alpha = alpha),
                SIMPLIFY = FALSE)
.lec.CurrentStreamEnd()

save(sims, file = paste0("./parts/design", dnum, "part", pnum, ".rda"))

