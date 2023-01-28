# Put simulation runs back together.
library(abind) # For concatenating arrays along one dimension

# Do this by design number, looping over the part numbers.
dnum <- as.numeric(commandArgs(trailingOnly = TRUE))

add <- function(l1, l2) {
  ans <- vector(mode = "list")
  ans$e_coef <- rbind(l1$e_coef, l2$e_coef)
  ans$e_se <- rbind(l1$e_se, l2$e_se)
  ans$e_cover <- rbind(l1$e_cover, l2$e_cover)

  ans$m_coef <- abind(l1$m_coef, l2$m_coef, along = 1)
  ans$m_t <- abind(l1$m_t, l2$m_t, along = 1)
  ans$m_mean <- rbind(l1$m_mean, l2$m_mean)
  ans$m_se <- rbind(l1$m_se, l2$m_se)
  ans$m_pCI <- abind(l1$m_pCI, l2$m_pCI, along = 1)
  ans$m_pCI2 <- abind(l1$m_pCI2, l2$m_pCI2, along = 1)
  ans$m_seCI <- abind(l1$m_seCI, l2$m_seCI, along = 1)
  ans$m_tCI <- abind(l1$m_tCI, l2$m_tCI, along = 1)
  ans$m_cover_p <- rbind(l1$m_cover_p, l2$m_cover_p)
  ans$m_cover_p2 <- rbind(l1$m_cover_p2, l2$m_cover_p2)
  ans$m_cover_se <- rbind(l1$m_cover_se, l2$m_cover_se)
  ans$m_cover_t <- rbind(l1$m_cover_t, l2$m_cover_t)

  ans$c_coef <- abind(l1$c_coef, l2$c_coef, along = 1)
  ans$c_t <- abind(l1$c_t, l2$c_t, along = 1)
  ans$c_mean <- rbind(l1$c_mean, l2$c_mean)
  ans$c_se <- rbind(l1$c_se, l2$c_se)
  ans$c_pCI <- abind(l1$c_pCI, l2$c_pCI, along = 1)
  ans$c_pCI2 <- abind(l1$c_pCI2, l2$c_pCI2, along = 1)
  ans$c_seCI <- abind(l1$c_seCI, l2$c_seCI, along = 1)
  ans$c_tCI <- abind(l1$c_tCI, l2$c_tCI, along = 1)
  ans$c_cover_p <- rbind(l1$c_cover_p, l2$c_cover_p)
  ans$c_cover_p2 <- rbind(l1$c_cover_p2, l2$c_cover_p2)
  ans$c_cover_se <- rbind(l1$c_cover_se, l2$c_cover_se)
  ans$c_cover_t <- rbind(l1$c_cover_t, l2$c_cover_t)

  ans$u_coef <- abind(l1$u_coef, l2$u_coef, along = 1)
  ans$u_t <- abind(l1$u_t, l2$u_t, along = 1)
  ans$u_mean <- rbind(l1$u_mean, l2$u_mean)
  ans$u_se <- rbind(l1$u_se, l2$u_se)
  ans$u_pCI <- abind(l1$u_pCI, l2$u_pCI, along = 1)
  ans$u_pCI2 <- abind(l1$u_pCI2, l2$u_pCI2, along = 1)
  ans$u_seCI <- abind(l1$u_seCI, l2$u_seCI, along = 1)
  ans$u_tCI <- abind(l1$u_tCI, l2$u_tCI, along = 1)
  ans$u_cover_p <- rbind(l1$u_cover_p, l2$u_cover_p)
  ans$u_cover_p2 <- rbind(l1$u_cover_p2, l2$u_cover_p2)
  ans$u_cover_se <- rbind(l1$u_cover_se, l2$u_cover_se)
  ans$u_cover_t <- rbind(l1$u_cover_t, l2$u_cover_t)

  ans$n_coef <- abind(l1$n_coef, l2$n_coef, along = 1)
  ans$n_t <- abind(l1$n_t, l2$n_t, along = 1)
  ans$n_mean <- rbind(l1$n_mean, l2$n_mean)
  ans$n_se <- rbind(l1$n_se, l2$n_se)
  ans$n_pCI <- abind(l1$n_pCI, l2$n_pCI, along = 1)
  ans$n_pCI2 <- abind(l1$n_pCI2, l2$n_pCI2, along = 1)
  ans$n_seCI <- abind(l1$n_seCI, l2$n_seCI, along = 1)
  ans$n_tCI <- abind(l1$n_tCI, l2$n_tCI, along = 1)
  ans$n_cover_p <- rbind(l1$n_cover_p, l2$n_cover_p)
  ans$n_cover_p2 <- rbind(l1$n_cover_p2, l2$n_cover_p2)
  ans$n_cover_se <- rbind(l1$n_cover_se, l2$n_cover_se)
  ans$n_cover_t <- rbind(l1$n_cover_t, l2$n_cover_t)

  ans$taus <- l1$taus # reports quantile levels
  ans$true_beta <- l1$true_beta # reports true coefficients
  ans$n <- l1$n
  ans$T <- l1$T
  ans$alpha <- l1$alpha
  ans$nml <- l1$nml
  ans$het <- l1$het
  ans$dep <- l1$dep
  ans$dyn <- l1$dyn
  return(ans)
}

load(paste0("./parts/design", dnum, "part1.rda")) # loads list sims
results <- sims

for (run in 2:100) {
  load(paste0("./parts/design", dnum, "part", run, ".rda"))
  results <- add(results, sims)
}
save(results, file = paste0("./design", dnum, "all.rda"))

