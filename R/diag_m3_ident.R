# Identification of Omega in m3 (stan/m3_corr.stan), from the full-sample fit
# output/fit_m3.rds (R/02_fit_models.R). Research only.
#
#   Rscript R/diag_m3_ident.R      # prints the checks, writes output/tables/m3_identification.csv
#
# With rho near 0, the data pin down v * 11' + Sigma (v = tau^2 / (1 - rho^2), the latent
# variance; Sigma = diag(sigma) Omega diag(sigma)) and the lag-1 autocovariance rho * v. The
# split between a common latent movement (v) and correlated measurement errors (Sigma_jk)
# is only weakly identified. The table records the evidence quoted in NOTES.md
# ("m3: rewrite and identification").

suppressPackageStartupMessages({ library(posterior) })
fit <- readRDS("output/fit_m3.rds")
d <- as_draws_df(fit$draws(c("c", "rho", "tau", "sigma", "Omega")))
v <- d$tau^2 / (1 - d$rho^2)
Sjk <- function(j, k) d[[sprintf("sigma[%d]", j)]] * d[[sprintf("sigma[%d]", k)]] * d[[sprintf("Omega[%d,%d]", j, k)]]
rows <- list()
add <- function(check, value) rows[[length(rows) + 1]] <<- data.frame(check = check, value = value)
add("v = tau^2/(1-rho^2): posterior mean", mean(v))
add("v: posterior sd", sd(v))
for (jk in list(c(1, 1), c(2, 2), c(3, 3), c(2, 1), c(3, 1), c(3, 2))) {
  s <- Sjk(jk[1], jk[2])
  add(sprintf("cor(v, Sigma[%d,%d])", jk[1], jk[2]), cor(v, s))
  add(sprintf("sd(v + Sigma[%d,%d]) / sd(Sigma[%d,%d])", jk[1], jk[2], jk[1], jk[2]), sd(v + s) / sd(s))
}
add("rho * v (lag-1 autocovariance): posterior mean", mean(d$rho * v))
add("rho * v: posterior sd", sd(d$rho * v))
avgO <- (d$`Omega[2,1]` + d$`Omega[3,1]` + d$`Omega[3,2]`) / 3
add("cor(tau, mean off-diagonal Omega)", cor(d$tau, avgO))
lo <- d$tau < quantile(d$tau, 0.1); hi <- d$tau > quantile(d$tau, 0.9)
for (nm in c("Omega[2,1]", "Omega[3,1]", "Omega[3,2]")) {
  add(sprintf("%s: mean in lowest-tau decile", nm), mean(d[[nm]][lo]))
  add(sprintf("%s: mean in highest-tau decile", nm), mean(d[[nm]][hi]))
}
add("tau: mean in lowest-tau decile", mean(d$tau[lo]))
add("tau: mean in highest-tau decile", mean(d$tau[hi]))
# LKJ(eta) on a K = 3 correlation matrix: (r + 1) / 2 ~ Beta(a, a), a = eta - 1 + K / 2
a <- 4 - 1 + 3 / 2
prior_sd <- 2 * sqrt(a * a / ((2 * a)^2 * (2 * a + 1)))
add("prior sd of an Omega off-diagonal under LKJ(4)", prior_sd)
for (nm in c("Omega[2,1]", "Omega[3,1]", "Omega[3,2]"))
  add(sprintf("%s: posterior sd / prior sd", nm), sd(d[[nm]]) / prior_sd)
g <- readRDS("data/processed/gdp_growth.rds"); y <- as.matrix(g[, c("g_E", "g_I", "g_P")])
vr <- apply(fit$draws("mu", format = "draws_matrix"), 1, var) / var(rowMeans(y))
add("variance ratio: 5% quantile", unname(quantile(vr, 0.05)))
add("variance ratio: median", median(vr))
add("variance ratio: 95% quantile", unname(quantile(vr, 0.95)))
tab <- do.call(rbind, rows)
write.csv(tab, "output/tables/m3_identification.csv", row.names = FALSE)
print(format(tab, digits = 3), row.names = FALSE)
