# Calculate confidence interval average coverage probabilities for several CI
# types.

library(Hmisc)

nvec <- c(20, 50, 100)
nT_prop <- c(2, 1, 1/2, 1/4, 1/10)
Nvec <- rep(nvec, each = length(nT_prop))
Tvec <- Nvec / nT_prop

taunames <- c("1/4", "1/2", "3/4")
method_names <- c("AT", "RWB.p", "RWB.se", "RWB.t")
methods_se <- c("AT", "RWB")

dnum <- 1

maketables <- function(dnum) {
  datname <- paste0("./data/design", dnum, "all.rda")
  load(datname) # loads an object called results
  nt_cols <- cbind(Nvec, Tvec)

  # Coverage probabilities
  cover_e_tab <- t(sapply(lapply(results, "[[", "e_cover"), colMeans))
  cover_m_p_tab <- t(sapply(lapply(results, "[[", "m_cover_p"), colMeans))
  cover_m_se_tab <- t(sapply(lapply(results, "[[", "m_cover_se"), colMeans))
  cover_m_t_tab <- t(sapply(lapply(results, "[[", "m_cover_t"), colMeans))
  cover_tab <- do.call(cbind, list(cover_e_tab, cover_m_p_tab, cover_m_se_tab,
                        cover_m_t_tab))
  cover_tab <- round(cover_tab * 100, 2)
  cover_tab <- cbind(nt_cols, cover_tab)
  colnames(cover_tab) <- c("n", "T", rep(taunames, (ncol(cover_tab) - 2) /
                            length(taunames)))

  # For naming tables
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
                length(method_names))), n.rgroup = rep(length(nT_prop),
                length(nvec)), caption = paste(desc, "Nominal 95 percent CIs."))
}

lapply(1:10, maketables)
