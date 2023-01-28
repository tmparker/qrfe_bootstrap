# Put simulation runs back together.
library(abind) # For concatenating arrays along one dimension

# Do this by design number, looping over the part numbers.
dnum <- as.numeric(commandArgs(trailingOnly = TRUE))[1]

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
  ans$m_seCI <- abind(l1$m_seCI, l2$m_seCI, along = 1)
  ans$m_tCI <- abind(l1$m_tCI, l2$m_tCI, along = 1)
  ans$m_cover_p <- rbind(l1$m_cover_p, l2$m_cover_p)
  ans$m_cover_se <- rbind(l1$m_cover_se, l2$m_cover_se)
  ans$m_cover_t <- rbind(l1$m_cover_t, l2$m_cover_t)
  ans$g_coef <- abind(l1$g_coef, l2$g_coef, along = 1)
  ans$g_t <- abind(l1$g_t, l2$g_t, along = 1)
  ans$g_mean <- rbind(l1$g_mean, l2$g_mean)
  ans$g_sd <- rbind(l1$g_sd, l2$g_sd)
  ans$g_pCI <- abind(l1$g_pCI, l2$g_pCI, along = 1)
  ans$g_seCI <- abind(l1$g_seCI, l2$g_seCI, along = 1)
  ans$g_tCI <- abind(l1$g_tCI, l2$g_tCI, along = 1)
  ans$g_cover_p <- rbind(l1$g_cover_p, l2$g_cover_p)
  ans$g_cover_se <- rbind(l1$g_cover_se, l2$g_cover_se)
  ans$g_cover_t <- rbind(l1$g_cover_t, l2$g_cover_t)
  ans$taus <- l1$taus
  ans$true_beta <- l1$true_beta
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

