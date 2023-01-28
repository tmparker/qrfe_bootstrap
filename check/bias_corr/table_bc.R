# Small table of coverage probabilities wiith and without bias correction

library(Hmisc)
load("design1all.rda") # loads an object called results

cover_e_tab <- colMeans(results$e_cover)
cover_s_tab <- colMeans(results$s_cover)
cover_sbc_tab <- colMeans(results$sbc_cover)
cover_j_tab <- colMeans(results$j_cover)
cover_m_p_tab <- colMeans(results$m_cover_p)
cover_m_se_tab <- colMeans(results$m_cover_se)
cover_m_t_tab <- colMeans(results$m_cover_t)
cover_tab <- cbind(cover_e_tab, cover_s_tab, cover_sbc_tab, cover_j_tab,
                    cover_m_p_tab, cover_m_se_tab, cover_m_t_tab)

# Heteroskedastic, dependent chi square error.

tau_names <- c("1/4", "1/2", "3/4")
method_names <- c("Standard estimator", "Smoothed estimator", "Smoothed,
  bias-corrected", "Smoothed, jackknife BC", "RWB.p", "RWB.se", "RWB.t")
mess <- "Empirical coverage of nominal 95 percent CIs.  The smoothed estimator
and bias corrections used correspond to those proposed in \\citet{GalvaoKato16}.
The RWB.x labels refer to the bootstrap confidence intervals considered in this
paper.  500 simulation runs with the heteroskedastic, dependent chi square
error design and \\(n = T = 100\\)."

cover_tab <- t(round(cover_tab * 100, 2))
dimnames(cover_tab) <- list(method_names, tau_names)
ctab <- Hmisc::latex(cover_tab, title = "", file = "bc_cover.tex", label = "bc_cover",
              size = "scriptsize", caption = mess)
