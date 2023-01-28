# Functions used by the smoothed estimators, also the bias correction estimator,
# which should only be used with the smoothed estimator but is applied to the
# standard estimator here.

# Kernel like Horowitz (1998)
K <- function(x) {
  k <- (105/64) * (1 - 5 * x^2 + 7 * x^4 - 3 * x^6)
  k[abs(x) > 1] <- 0
  k
}

dK <- function(x) {
  dk <- (105/64) * (-10 * x + 28 * x^3 - 18 * x^5)
  dk[abs(x) > 1] <- 0
  dk
}

# Galvao & Kato (2016) integrated kernel
G <- function(x) {
  g <- 0.5 - (105/64) * (x - 5/3 * x^3 + 7/5 * x^5 - 3/7 * x^7)
  g[x < -1] <- 1
  g[x > 1] <- 0
  g
}

# Smoothed QR objective function
# Assumes x is a matrix
srq_opt <- function(coefs, tau, y, x, h) {
  uhat <- y - x %*% coefs
  sum(uhat * (tau - G(uhat / h)))
}

# Smoothed (obj. fcn.) coefficient estimator
srq <- function(y, x, d, rq_obj, n, T, tau) {
  h <- sd(rq_obj$resid) * (n * T)^(-1/7) # G&K p.99
  fit <- optim(rq_obj$coef, srq_opt, tau = tau,
                y = y, x = cbind(x, d), h = h)
  coefs <- fit$par
  uhat <- y - coefs[1] * x - rep(coefs[2:(n+1)], each = T)
  list(coefficients = coefs, residuals = uhat)
}

# Conquer (= smoothed estimating equation) estimator
conq <- function(y, x, d, n, T, tau) {
  dm1 <- d[, -n] # conquer() has to include an intercept
  fit <- conquer(cbind(x, dm1), y, tau = tau, ci = "none")
  coefs <- fit$coef[2]
  uhat <- fit$resid
  list(coefficients = coefs, residuals = uhat)
}

# Bias correction for smoothed estimator (also applied to standard estimator)
srq_biascorrect <- function(y, x, d, n, T, tau) {
  standard <- rq.fit(cbind(x, d), y, tau)
  betahat <- standard$coef[1]
  smoothed <- srq(y, x, d, standard, n, T, tau)
  uhat <- smoothed$resid
  bn <- 2 * sd(uhat) * T^(-1/5)
  mn <- 1 # G&K 2016 choice
  # Galvao & Kato p.97
  ind <- rep(1:n, each = T)
  fi <- tapply(K(uhat / bn), ind, mean) / bn
  si <- double(n)
  si[fi > cc] <- fi^(-1)
  gammai <- tapply(K(uhat / bn) * x, ind, mean) * si / bn
  xmg <- x - rep(gammai, each = T)
  nui <- tapply(dK(uhat / bn) * xmg, ind, mean) / bn^2
  Gam <- mean(K(uhat / bn) * x * xmg) / bn
  # For omegas, each unit has an estimate from -mn to +mn
  ph <- vph <- vrh <- matrix(0, n, 2 * mn + 1)
  i1 <- c((mn + 1):1, rep(1, mn)) # to index leads/lags
  i2 <- c(rep(T, mn), T:(T - mn))
  i3 <- rev(i1)
  i4 <- rev(i2)
  for (i in 1:n) {
    for (j in 1:(2 * mn + 1)) {
      ku.ind <- (i - 1) * T + i1[j]:i2[j]
      ulow.ind <- (i - 1) * T + i3[j]:i4[j]
      ph[i, j] <- sum( K(uhat / bn)[ku.ind] * (uhat <= 0)[ulow.ind] ) / T
      vph[i, j] <- sum( K(uhat / bn)[ku.ind] * x[ku.ind] * (uhat <= 0)[ulow.ind] ) / T
      vrh[i, j] <- sum( (uhat <= 0)[ku.ind] * (uhat <= 0)[ulow.ind] ) / T
    }
  }
  # Their omegas are weighted averages of the above things
  w_seq <- 1 - abs((-mn):mn) / (mn + 1)
  w_seq[mn + 1] <- 0
  om1 <- rowSums(w_seq * sweep(-ph, 1, tau * fi, FUN = "+"))
  om2 <- rowSums(w_seq * sweep(-vph, 1, tau * fi * gammai, FUN = "+"))
  om3 <- tau * (1 - tau) + rowSums(w_seq * sweep(vrh, 1, tau^2))
  bhat <- mean(si * (om1 * gammai - om2 + 0.5 * si * om3 * nui)) / Gam
  # bsqr <- betahat_smooth - bhat / T
  bsqr <- smoothed$coef[1] - bhat / T
  bqr <- betahat - bhat / T
  ests <- c(bqr, bsqr)
  names(ests) <- c("bcrq", "bcsrq")
  ests
}

half_panel_jackknife <- function(y, x, d, n, T, tau) {
  standard <- rq.fit(cbind(x, d), y, tau)
  betahat_smooth <- srq(y, x, d, standard, n, T, tau)
  T2c <- ceiling(T / 2)
  T2f <- floor(T / 2)
  half_ind <- rep(c(rep(1, T2c), rep(2, T2f)), n)
  y11 <- y[half_ind == 1]
  y12 <- y[half_ind == 2]
  x11 <- x[half_ind == 1]
  x12 <- x[half_ind == 2]
  d11 <- kronecker(diag(n), rep(1, T2c))
  d12 <- kronecker(diag(n), rep(1, T2f))
  start11 <- rq.fit(cbind(x11, d11), y11, tau)
  start12 <- rq.fit(cbind(x12, d12), y12, tau)
  hn <- sd(start11$resid) * (n * T)^(-1/7)
  fit11 <- srq(y11, x11, d11, start11, n = n, T = T2c, tau = tau)
  fit12 <- srq(y12, x12, d12, start12, n = n, T = T2f, tau = tau)
  beta11 <- fit11$coef[1]
  beta12 <- fit12$coef[1]
  bbar1 <- (T2c / T) * beta11 + (T2f / T) * beta12
  bbar <- bbar1
  if (T %% 2 != 0) {
    half_ind2 <- rep(c(rep(1, T2f), rep(2, T2c)), n)
    y21 <- y[half_ind2 == 1]
    y22 <- y[half_ind2 == 2]
    x21 <- x[half_ind2 == 1]
    x22 <- x[half_ind2 == 2]
    start21 <- rq.fit(cbind(x21, d12), y21, tau) # reverse ind.
    start22 <- rq.fit(cbind(x22, d11), y22, tau)
    fit21 <- srq(y21, x21, d12, start21, n = n, T = T2f, tau = tau)
    fit22 <- srq(y22, x22, d11, start22, n = n, T = T2c, tau = tau)
    beta21 <- fit21$coef[1]
    beta22 <- fit22$coef[1]
    bbar2 <- (T2f / T) * beta21 + (T2c / T) * beta22
    bbar <- (bbar1 + bbar2) / 2
  }
  2 * betahat_smooth$coef[1] - bbar
}

