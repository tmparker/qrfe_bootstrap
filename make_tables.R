# Take a set of results and make some tables by rearranging the entries.

library(Hmisc)

# dnum <- as.numeric(commandArgs(trailingOnly = TRUE))[1]
#dnum <- 1

nvec <- c(25, 50, 100)
tvec <- c(20, 50, 100)

length_expand <- 100
bias_expand <- 1e5
rmse_expand <- 1e5

taunames <- c("1/4", "1/2", "3/4")
# method_names <- c("Asy.", "M.p", "M.se", "M.t", "G.p", "G.se", "G.t")
method_names <- c("AT", "RWB.p", "RWB.se", "RWB.t")
#methods_se <- c("Asy.", "Mult.", "Grad.")
#methods_se <- c("Asy.", "Mult.")
methods_se <- c("AT", "RWB")
nt_cols <- cbind(rep(nvec, each = length(tvec)), rep(tvec, length(nvec)))

maketables <- function(dnum) {
datname <- paste0("./data/design", dnum, "all.rda")
load(datname) # loads an object called results

# Coverage probabilities
cover_e_tab <- t(sapply(lapply(results, "[[", "e_cover"), colMeans))
cover_m_p_tab <- t(sapply(lapply(results, "[[", "m_cover_p"), colMeans))
cover_m_se_tab <- t(sapply(lapply(results, "[[", "m_cover_se"), colMeans))
cover_m_t_tab <- t(sapply(lapply(results, "[[", "m_cover_t"), colMeans))
# cover_g_p_tab <- t(sapply(lapply(results, "[[", "g_cover_p"), colMeans))
# cover_g_se_tab <- t(sapply(lapply(results, "[[", "g_cover_se"), colMeans))
# cover_g_t_tab <- t(sapply(lapply(results, "[[", "g_cover_t"), colMeans))
# cover_tab <- do.call(cbind, list(cover_e_tab, cover_m_p_tab, cover_m_se_tab,
#                       cover_m_t_tab, cover_g_p_tab, cover_g_se_tab, cover_g_t_tab))
cover_tab <- do.call(cbind, list(cover_e_tab, cover_m_p_tab, cover_m_se_tab,
                      cover_m_t_tab))
cover_tab <- round(cover_tab * 100, 2)
cover_tab <- cbind(nt_cols, cover_tab)
colnames(cover_tab) <- c("n", "T", rep(taunames, (ncol(cover_tab) - 2) /
                          length(taunames)))

# Average CI length
z_al2 <- qnorm(1 - results[[1]]$alpha / 2)
length_e <- lapply(lapply(results, function(x) x$e_sd), function(y) y * z_al2 * 2)
length_e_tab <- t(sapply(length_e, colMeans))
length_m_p_tab <- t(sapply(results, function(x) rowMeans(apply(x$m_pCI, 1, diff))))
length_m_se_tab <- t(sapply(results, function(x) rowMeans(apply(x$m_seCI, 1, diff))))
length_m_t_tab <- t(sapply(results, function(x) rowMeans(apply(x$m_tCI, 1, diff))))
# length_g_p_tab <- t(sapply(results, function(x) rowMeans(apply(x$g_pCI, 1, diff))))
# length_g_se_tab <- t(sapply(results, function(x) rowMeans(apply(x$g_seCI, 1, diff))))
# length_g_t_tab <- t(sapply(results, function(x) rowMeans(apply(x$g_tCI, 1, diff))))
# length_tab <- do.call(cbind, list(length_e_tab, length_m_p_tab,
#                         length_m_se_tab, length_m_t_tab, length_g_p_tab,
#                         length_g_se_tab, length_g_t_tab))
length_tab <- do.call(cbind, list(length_e_tab, length_m_p_tab,
                        length_m_se_tab, length_m_t_tab))
length_tab <- round(length_tab * length_expand, 1)
length_tab <- cbind(nt_cols, length_tab)
colnames(length_tab) <- c("n", "T", rep(taunames, (ncol(length_tab) - 2) /
                          length(taunames)))

# Bias of variance estimates
v <- t(sapply(results, "[[", "true_se"))^2
ev <- t(sapply(results, function(x) colMeans(x$e_sd^2)))
mv <- t(sapply(results, function(x) colMeans(x$m_sd^2)))
# gv <- t(sapply(results, function(x) colMeans(x$g_sd^2)))
# v_bias_tab <- do.call(cbind, list(ev - v, mv - v, gv - v))
v_bias_tab <- do.call(cbind, list(ev - v, mv - v))
v_bias_tab <- round(v_bias_tab * bias_expand, 1)
v_bias_tab <- cbind(nt_cols, v_bias_tab)
colnames(v_bias_tab) <- c("n", "T", rep(taunames, (ncol(v_bias_tab) - 2) /
                          length(taunames)))

# RMSE of variance estimates
e_v_rmse <- m_v_rmse <- g_v_rmse <- matrix(0, 9, 3)
for (i in 1:9) { # length of nvec times length of tvec
  e_v_rmse[i, ] <- sqrt(rowMeans((t(results[[i]]$e_sd^2) - v[i, ])^2))
  m_v_rmse[i, ] <- sqrt(rowMeans((t(results[[i]]$m_sd^2) - v[i, ])^2))
  # g_v_rmse[i, ] <- sqrt(rowMeans((t(results[[i]]$g_sd^2) - v[i, ])^2))
}
# v_rmse_tab <- do.call(cbind, list(e_v_rmse, m_v_rmse, g_v_rmse))
v_rmse_tab <- do.call(cbind, list(e_v_rmse, m_v_rmse))
v_rmse_tab <- round(v_rmse_tab * rmse_expand, 1)
v_rmse_tab <- cbind(nt_cols, v_rmse_tab)
colnames(v_rmse_tab) <- c("n", "T", rep(taunames, (ncol(v_rmse_tab) - 2) /
                          length(taunames)))

# For naming tables - unfortunately later it was decided that
# heteroskedasticity/homoskedasticity should come before dependent/independent
# in the naming convention.
if (dnum == 1) {file_prefix <- "chisq_dep_het"
} else if (dnum == 2) {file_prefix <- "chisq_dep_loc"
} else if (dnum == 3) {file_prefix <- "chisq_dyn"
} else if (dnum == 4) {file_prefix <- "chisq_ind_het"
} else if (dnum == 5) {file_prefix <- "chisq_ind_loc"
} else if (dnum == 6) {file_prefix <- "norm_dep_het"
} else if (dnum == 7) {file_prefix <- "norm_dep_loc"
} else if (dnum == 8) {file_prefix <- "norm_dyn"
} else if (dnum == 9) {file_prefix <- "norm_ind_het"
} else {file_prefix <- "norm_ind_loc"}

# For describing DGPs in table captions
if (dnum == 1) {desc <- "Heteroskedastic, dependent chi square error."
} else if (dnum == 2) {desc <- "Homoskedastic, dependent chi square error."
} else if (dnum == 3) {desc <- "Dynamic chi square model."
} else if (dnum == 4) {desc <- "Heteroskedastic, independent chi square error."
} else if (dnum == 5) {desc <- "Homoskedastic, independent chi square error."
} else if (dnum == 6) {desc <- "Heteroskedastic, dependent normal error."
} else if (dnum == 7) {desc <- "Homoskedastic, dependent normal error."
} else if (dnum == 8) {desc <- "Dynamic normal model."
} else if (dnum == 9) {desc <- "Heteroskedastic, independent normal error."
} else {desc <- "Homoskedastic, independent normal error."}

ctab <- latex(cover_tab, title = "",
              file = paste0("./tables/", file_prefix, "_cover.tex"),
              label = paste0(file_prefix, "_cover"), size = "scriptsize",
              cgroup = c("", method_names), n.cgroup = c(2, rep(length(taunames),
              (ncol(cover_tab) - 2) / length(taunames))),
              n.rgroup = rep(length(tvec), length(nvec)),
              caption = paste(desc, "Nominal 90 percent CIs, multiplier
                bootstrap."))

ltab <- latex(length_tab, title = "", file = paste0("./tables/",
      file_prefix, "_length.tex"), label = paste0(file_prefix, "_length"),
      size = "scriptsize", cgroup = c("", method_names),
      n.cgroup = c(2, rep(length(taunames), (ncol(length_tab) - 2) /
      length(taunames))), n.rgroup = rep(length(tvec), length(nvec)),
      caption = paste(desc, "Average CI lengths x 100, multiplier bootstrap."))

btab <- latex(v_bias_tab, title = "Bias of var ests.", file =
      paste0("./tables/", file_prefix, "_bias.tex"), label = paste0(file_prefix,
      "_bias"), size = "scriptsize", cgroup = c("", methods_se),
      n.cgroup = c(2, rep(length(taunames), (ncol(v_bias_tab) - 2) /
      length(taunames))), n.rgroup = rep(length(tvec), length(nvec)),
      caption = paste(desc, "Bias of variance estimates x 1e5."))

mtab <- latex(v_rmse_tab, title = "RMSE of var ests.", file =
      paste0("./tables/", file_prefix, "_rmse.tex"), label = paste0(file_prefix,
      "_rmse"), size = "scriptsize", cgroup = c("", methods_se),
      n.cgroup = c(2, rep(length(taunames), (ncol(v_rmse_tab) - 2) /
      length(taunames))), n.rgroup = rep(length(tvec), length(nvec)),
      caption = paste(desc, "RMSE of variance estimates x 1e5."))
}

lapply(1:10, maketables)
