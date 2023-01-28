# These functions and parameters are used in all simulations

# These parameters are constant across all the simulations
taus <- 1:3 / 4
al <- 0.05
simreps <- 5

cc <- 0.01 #Trimming for sparsity estimation
nboot <- 999
nvec <- c(20, 50, 100)

# List of designs
# Designs are numbered 1-10, baed on alphabetical order of names.  For example,
# chisq_dep_het is the first one alphabetically, so it's design 1.  The result
# of this function is used to set the arguments for the subsequent functions.
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
dgp <- function(n, T, nml = TRUE, het = FALSE, dep = FALSE, dyn = FALSE,
                be = 1, ga = 0.2, the = 0.5, rho = 0.4, dyn_burnin = 50) {
  if (dyn) {
    al <- runif(n)
    u <- matrix(rnorm(n * (T + dyn_burnin)), T + dyn_burnin, n)
    yaux <- matrix(0, T + dyn_burnin, n)
    yaux[1, ] <- al + u[1, ]
    for(i in 2:(T + dyn_burnin)) {
      yaux[i, ] <- al + rho * yaux[i-1, ] + u[i, ]
    }
    yaux.y <- yaux[(dyn_burnin + 1):(T + dyn_burnin), ]
    yaux.yl <- yaux[dyn_burnin:(T + dyn_burnin-1), ]
    y <- c(yaux.y)
    x <- c(yaux.yl)
    return(list(y = y, x = x, true_al = al))
  } else {
    if (!het) { ga <- 0 }
    al <- rep(runif(n), each = T)
    x <- rchisq(n * T, df = 3) + 0.3 * al
    if (dep) { # MA(1) error terms if dependent:
      if (nml) {
        e <- rnorm(n * (T + 1))
      } else {
        e <- rchisq(n * (T + 1), df = 4)
      }
      i1 <- rep((1:n - 1) * (T + 1), each = T) + 2:(T + 1)
      i0 <- rep((1:n - 1) * (T + 1), each = T) + 1:T
      u <- e[i1] + the * e[i0]
    } else { # for independent errors
      if (nml) {
        u <- rnorm(n * T)
      } else {
        u <- rchisq(n * T, df = 4)
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

# Asymptotic standard error estimator using the defaults described in Kato,
# Galvao & Montes-Rojas (2012) and Galvao & Kato (2016)
est_se <- function(coefs, uhat, x, n, T, tau, indep = FALSE) {
  x <- as.matrix(x)
  p <- ncol(x)
  coefs <- coefs[1:p]
  ind <- rep(1:n, each = T)
  m <- 1 # bandwidth for lag lengths
  scl <- min(sd(uhat), (quantile(uhat, 0.75) - quantile(uhat, 0.25)) / 1.34)
  h <- 2 * scl * T^-0.2
  f <- dnorm(uhat / h) / h
  fi <- tapply(f, ind, mean)
  si <- double(n)
  si[fi > cc] <- fi[fi > cc]^(-1)
  gammai <- apply(f * as.matrix(x), 2, function(M) tapply(M, ind, mean)) * si
  xmg <- x - rep(gammai, each = T)
  Ginv <- diag(p)
  Ginv <- backsolve(qr(sqrt(f) * xmg)$qr[1:p, 1:p, drop = FALSE], Ginv)
  Ginv <- tcrossprod(Ginv)
  if (indep == TRUE) {
    V <- tau * (1 - tau) * crossprod(xmg)
  } else {
    w_seq <- 1 - abs((-m):m) / (m + 1)
    w_seq[m + 1] <- 0
    i1 <- c((m + 1):1, rep(1, m)) # to index leads/lags, could use symmetry
    i2 <- c(rep(T, m), T:(T - m))
    i3 <- rev(i1)
    i4 <- rev(i2)
    vi_list <- vector(mode = "list", length = n)
    for (i in 1:n) {
      vi_list[[i]] <- vector(mode = "list", length = 2 * m + 1)
      for (j in 1:(2 * m + 1)) {
        ind_now <- (i - 1) * T + i1[j]:i2[j]
        ind_lag <- (i - 1) * T + i3[j]:i4[j]
        part_now <- (tau - (uhat <= 0)[ind_now]) * xmg[ind_now, , drop = FALSE]
        part_lag <- (tau - (uhat <= 0)[ind_lag]) * xmg[ind_lag, , drop = FALSE]
        vi_list[[i]][[j]] <- w_seq[j] * crossprod(part_now, part_lag)
      }
      vi_list[[i]][[m + 1]] <- tau * (1 - tau) *
                                crossprod(xmg[(i - 1) * T + 1:T, , drop = FALSE])
    }
    Vi <- lapply(vi_list, function(L) Reduce('+', L))
    V <- Reduce('+', Vi)
  }
  GVG <- Ginv %*% V %*% Ginv
  sqrt(diag(GVG))
}

