# Compare empirical CI coverage probabilities for two percentile CI methods

library(Hmisc)
load("design1all.rda") # loads an object called results

cover_tab <- cbind(colMeans(results$e_cover), colMeans(results$m_cover_p1),
                    colMeans(results$m_cover_p2))
cover_tab <- round(cover_tab * 100, 2)

taunames <- c("1/4", "1/2", "3/4")
method_names <- c("Standard estimator", "RWB percentile 1", "RWB percentile 2")
dimnames(cover_tab) <- list(method_names, taunames)

desc <- "Heteroskedastic, dependent chi square error."
mess <- "Empirical coverage of nominal 95 percent CIs.  Comparison of two
bootstrap percentile CIs.  Letting \\(\\widehat{\\beta}^\\ast_\\eta\\) be the
\\(\\eta\\)-th percentile of the bootstrap distribution of
\\(\\widehat{\\beta}^\\ast\\), method 1 is a CI of the form \\(
  (\\widehat{\\beta} - \\widehat{\\beta}^\\ast_{1 - \\alpha / 2},
  \\widehat{\\beta} - \\widehat{\\beta}^\\ast_{\\alpha / 2}) \\), while method
  2 is a CI of the form \\( (\\widehat{\\beta}^\\ast_{\\alpha / 2},
  \\widehat{\\beta}^\\ast_{1 - \\alpha / 2}) \\).  1000 simulation runs with
  the heteroskedastic, dependent chi square error design and \\(n = T = 100\\).
  Method 2 was used in the main simulations."

ctab <- latex(cover_tab, title = "", file = "perc_cover.tex", label =
  "perc_cover", size = "scriptsize", caption = mess)

