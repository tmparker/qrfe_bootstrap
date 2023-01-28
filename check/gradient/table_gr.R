# Coverage probabilities for objective-multiplier vs. gradient-multiplier
# bootstraps.

library(Hmisc)
load("design1all.rda") # loads an object called results

# Coverage probabilities
cover_e_tab <- colMeans(results$e_cover)
cover_m_p_tab <- colMeans(results$m_cover_p)
cover_m_se_tab <- colMeans(results$m_cover_se)
cover_m_t_tab <- colMeans(results$m_cover_t)
cover_g_p_tab <- colMeans(results$g_cover_p)
cover_g_se_tab <- colMeans(results$g_cover_se)
cover_g_t_tab <- colMeans(results$g_cover_t)
cover_tab <- cbind(cover_e_tab, cover_m_p_tab, cover_m_se_tab,
                      cover_m_t_tab, cover_g_p_tab, cover_g_se_tab, cover_g_t_tab)
cover_tab <- t(round(cover_tab * 100, 2))

tau_names <- c("1/4", "1/2", "3/4")
method_names <- c("AT", "RWB.p", "RWB.se", "RWB.t", "G.p", "G.se", "G.t")
mess <- "Nominal 95 percent CIs.  Objective function multiplier and gradient
multiplier bootstraps.  AT stands for the standard estimator using
asymptotic theory, while RWB.x and G.x stand for random weighted bootstrap and
gradient bootstrap CIs, where the gradient bootstrap is from
\\citet{Hagemann17}.  1000 simulation runs with the heteroskedastic, dependent
chi square error design and \\(n = T = 100\\)."
dimnames(cover_tab) <- list(method_names, tau_names)

ctab <- Hmisc::latex(cover_tab, title = "", file = "grad_cover.tex",
              label = "grad_cover", size = "scriptsize", caption = mess)
