# Put simulation runs back together.
# library(abind) # For concatenating arrays along one dimension

# Do this by design number, looping over the part numbers.
dnum <- as.numeric(commandArgs(trailingOnly = TRUE))

add <- function(l1, l2) {
  ans <- vector(mode = "list")
  ans$e_coef <- rbind(l1$e_coef, l2$e_coef)
  ans$s_coef <- rbind(l1$s_coef, l2$s_coef)
  ans$c_coef <- rbind(l1$c_coef, l2$c_coef)
  ans$e_se <- rbind(l1$e_se, l2$e_se)
  ans$s_se <- rbind(l1$s_se, l2$s_se)
  ans$c_se <- rbind(l1$c_se, l2$c_se)
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

