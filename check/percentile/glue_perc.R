# Put simulation runs together.
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
  ans$m_pCI1 <- abind(l1$m_pCI1, l2$m_pCI1, along = 1)
  ans$m_pCI2 <- abind(l1$m_pCI2, l2$m_pCI2, along = 1)
  ans$m_cover_p1 <- rbind(l1$m_cover_p1, l2$m_cover_p1)
  ans$m_cover_p2 <- rbind(l1$m_cover_p2, l2$m_cover_p2)
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

