# Looking at some features of different bootstrap weight distributions

load("design1all.rda") # loads list results
taus <- 1:3 / 4
ntau <- length(taus)
namevec <- c("True", "exp.", "chi sq.", "unif.", "f. norm.")

# Boxplots of parameter estimate distributions
est_sel <- results[c("e_coef", "m_coef", "c_coef", "u_coef", "n_coef")]
est <- results$e_coef

source("../../common.R") # to get the function beta_true()
bt <- beta_true(taus, nml = FALSE, het = TRUE, dep = TRUE)
est_err <- t(t(est) - bt)

cw_df <- 5 # degrees of freedom parameter for chi square weights
sc <- sqrt(2 / cw_df) # SD of scaled chi square weights
su <- sqrt(1 / 3) # SD of U[0, 2] weights
sn <- sqrt(pi / 2 - 1) # SD of folded normal weights (s.t. E[W] = 1)
sc_vec <- c(1, 1, sc, su, sn) # se ests scaled in sim_weight.R

# Consider the first simulation repetition as representative
i <- 1
e_err <- t(t(est_sel$m_coef[i, , ]) - est[i, ])
c_err <- t(t(est_sel$c_coef[i, , ]) - est[i, ]) / sc
u_err <- t(t(est_sel$u_coef[i, , ]) - est[i, ]) / su
n_err <- t(t(est_sel$n_coef[i, , ]) - est[i, ]) / sn
err <- list(est_err, e_err, c_err, u_err, n_err)
names(err) <- namevec

pdf(file = "weight_est_boxplots.pdf", height = 5, width = 12)
par(mfrow = c(1, 3), mar = c(3.1, 3.1, 1.1, 1.1), oma = c(0, 0, 4, 0), mgp = 2:0)
for (k in 1:ntau) {
  boxplot(err, main = paste("tau =", taus[k]), notch = TRUE)
  abline(h = 0, col = "red", lty = 2)
}
mtext("Effect of bootstrap weight distributions on bootstrap parameter estimates",
      outer = TRUE, cex = 1.5)
dev.off()

# Boxplots of t distributions
est_t <- est_err / results$e_se
# Once again, use i=1 as representative
t_sel <- results[c("m_t", "c_t", "u_t", "n_t")]
t_sel <- lapply(t_sel, function(M) as.matrix(M[1, , ]))
t_for_box <- c(list(est_t = est_t), t_sel)
names(t_for_box) <- namevec
t_box_rescale <- mapply(function(M, s) M / s, t_for_box, sc_vec)

pdf(file = "weight_t_boxplots.pdf", height = 5, width = 12)
par(mfrow = c(1, 3), mar = c(3.1, 3.1, 1.1, 1.1), oma = c(0, 0, 4, 0), mgp = 2:0)
for (k in 1:ntau) {
  boxplot(t_box_rescale, main = paste("tau =", taus[k]), notch = TRUE)
  abline(h = 0, col = "red", lty = 2)
}
mtext("Effect of bootstrap weight distributions on t statistics",
      outer = TRUE, cex = 1.5)
dev.off()

# Boxplots of variance estimates
true_se <- apply(results$e_coef, 2, sd)
se_sel <- results[c("e_se", "m_se", "c_se", "u_se", "n_se")]
se_scl <- mapply(function(M, s) M / s, se_sel, sc_vec, SIMPLIFY = FALSE)
var_errors_scl <- lapply(se_scl, function(M) t(t(M^2) - true_se^2))
names(var_errors_scl) <- c("AT", "exp.", "chi sq.", "unif.", "f. norm.")

pdf(file = "weight_var_boxplots.pdf", height = 5, width = 12)
par(mfrow = c(1, 3), mar = c(3.1, 3.1, 1.1, 1.1), oma = c(0, 0, 4, 0), mgp = 2:0)
for (k in 1:ntau) {
  vk <- sapply(var_errors_scl, function(A) A[, k])
  boxplot(vk, main = paste("tau =", taus[k]), notch = TRUE)
  abline(h = 0, col = "red", lty = 2)
}
mtext("Effect of bootstrap weight distributions on variance estimation",
      outer = TRUE, cex = 1.5)
dev.off()

# Table of empirical coverage probabilities

method_names <- c("AT", "Exponential", "Chi square", "Uniform", "Folded normal")
tau_names <- c("1/4", "1/2", "3/4")
cover_tab_p <- cbind(colMeans(results$e_cover), colMeans(results$m_cover_p),
                      colMeans(results$c_cover_p), colMeans(results$u_cover_p),
                      colMeans(results$n_cover_p))
dimnames(cover_tab_p) <- list(tau_names, method_names)
cat("percentile (1) \n")
print(round(cover_tab_p * 100, 2))

cover_tab_p2 <- cbind(colMeans(results$e_cover), colMeans(results$m_cover_p2),
                      colMeans(results$c_cover_p2), colMeans(results$u_cover_p2),
                      colMeans(results$n_cover_p2))
dimnames(cover_tab_p2) <- list(tau_names, method_names)
cat("percentile (2) \n")
print(round(cover_tab_p2 * 100, 2))

cover_tab_se <- cbind(colMeans(results$e_cover), colMeans(results$m_cover_se),
                      colMeans(results$c_cover_se), colMeans(results$u_cover_se),
                      colMeans(results$n_cover_se))
dimnames(cover_tab_se) <- list(tau_names, method_names)
cat("standard error \n")
print(round(cover_tab_se * 100, 2))

cover_tab_t <- cbind(colMeans(results$e_cover), colMeans(results$m_cover_t),
                      colMeans(results$c_cover_t), colMeans(results$u_cover_t),
                      colMeans(results$n_cover_t))
dimnames(cover_tab_t) <- list(tau_names, method_names)
cat("t \n")
print(round(cover_tab_t * 100, 2))

# mess <- "Empirical coverage of nominal 95 percent CIs.  The smoothed estimator
# and bias corrections used correspond to those proposed in \\citet{GalvaoKato16}.
# The RWB.x labels refer to the bootstrap confidence intervals considered in this
# paper.  1000 simulation runs with the heteroskedastic, dependent chi square
# error design and \\(n = T = 100\\)."

# cover_tab_p <- t(round(cover_tab_p * 100, 2))
# dimnames(cover_tab_p) <- list(method_names, tau_names)
# ctab_p <- Hmisc::latex(cover_tab_p, title = "", file = "bc_cover.tex", label = "bc_cover",
#               size = "scriptsize", caption = mess)


