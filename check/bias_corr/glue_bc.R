# Put simulation runs back together.
library(abind) # For concatenating arrays along one dimension

# Do this by design number, looping over the part numbers.
dnum <- as.numeric(commandArgs(trailingOnly = TRUE))

add <- function(l1, l2) {
  ans <- vector(mode = "list")
  ans$e_coef <- rbind(l1$e_coef, l2$e_coef)
  ans$e_se <- rbind(l1$e_se, l2$e_se)
  ans$e_cover <- rbind(l1$e_cover, l2$e_cover)

  ans$s_coef <- rbind(l1$s_coef, l2$s_coef)
  ans$s_se <- rbind(l1$s_se, l2$s_se)
  ans$s_cover <- rbind(l1$s_cover, l2$s_cover)

  ans$ebc_coef <- rbind(l1$ebc_coef, l2$ebc_coef)
  ans$ebc_cover <- rbind(l1$ebc_cover, l2$ebc_cover)
  ans$sbc_coef <- rbind(l1$sbc_coef, l2$sbc_coef)
  ans$sbc_cover <- rbind(l1$sbc_cover, l2$sbc_cover)
  ans$j_coef <- rbind(l1$j_coef, l2$j_coef)
  ans$j_cover <- rbind(l1$j_cover, l2$j_cover)

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

