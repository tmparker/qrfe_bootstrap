# Side-by-side density plots of coefficient estimates, standard error estimates
# and t distributions with rq, srq and conquer estimators

load("design1all.rda")

# calculate t statistics
e_t <- t(t(results$e_coef) - results$true_beta) / results$e_se
s_t <- t(t(results$s_coef) - results$true_beta) / results$s_se
c_t <- t(t(results$c_coef) - results$true_beta) / results$c_se

pdf(file = "smooth_all.pdf", height = 8, width = 8)
par(mfrow = c(3, 3), mar = c(3.1, 3.1, 1.1, 1.1), mgp = 2:0)
for (k in 1:3) {
  plot(density(results$e_coef[, k]), main = "Coefficients")
  lines(density(results$s_coef[, k]), col = "red", lty = 2)
  lines(density(results$c_coef[, k]), col = "blue", lty = 3)
  abline(h = 0, lty = 2, col = "gray50")

  plot(density(results$e_se[, k]), main = "Standard err.")
  lines(density(results$s_se[, k]), col = "red", lty = 2)
  lines(density(results$c_se[, k]), col = "blue", lty = 3)
  abline(h = 0, lty = 2, col = "gray50")

  plot(density(e_t[, k]), main = "t stat.")
  lines(density(s_t[, k]), col = "red", lty = 2)
  lines(density(c_t[, k]), col = "blue", lty = 3)
  abline(h = 0, lty = 2, col = "gray50")

}
dev.off()


