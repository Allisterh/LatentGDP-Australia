# Research-suite figures from the saved fits -> <root>/figures/*.png.
#   Rscript R/04_figures.R
#   LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980 Rscript R/04_figures.R
# Sample window and output root: R/sample_window.R (the fits under <root> must have been
# made with the same LATENTGDP_START). Both unset = the full sample, output/, as before.
# Figure 07 is drawn only if <root>/tables/lfo_pointwise.csv exists (R/03_compare.R).

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

theme_set(theme_minimal(base_size = 11))

source("R/sample_window.R")
cat(sw_describe(), "\n")
fig_dir <- sw_out("figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
fit_path <- function(nm) sw_out(sprintf("fit_%s.rds", nm))
# Appended to the subtitles of figures 01, 02, 04, 06 and 07 when the sample is truncated
# (empty otherwise, so the full-sample figures are unchanged). Figure 05's subtitle already
# gives its sample start; figure 03 (2018 onward) has no subtitle.
samp_note <- if (is.null(sw_start())) "" else
  sprintf(" Sample from %s (LATENTGDP_START).", format(sw_start(), "%Y-%m"))

gdp_growth   <- sw_trim_q(readRDS("data/processed/gdp_growth.rds")) |>
  mutate(headline = (g_E + g_I + g_P) / 3)
gdp_growth_4 <- sw_trim_q(readRDS("data/processed/gdp_growth_with_hours.rds")) |>
  mutate(headline = (g_E + g_I + g_P) / 3)

names_full <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m9", "m12", "m13")
labels_m <- c(
  m1 = "m1: baseline", m2 = "m2: AR(1) latent",
  m3 = "m3: free corr errors", m4 = "m4: SV",
  m5 = "m5: serial errors", m6 = "m6: structural Sigma",
  m7 = "m7: COVID scale", m8 = "m8: 4-signal (+hours)",
  m9 = "m9: structural Sigma + SV", m10 = "m10: 4-signal + SV (latent SV)",
  m11 = "m11: mixed-frequency", m12 = "m12: m9 + serial errors",
  m13 = "m13: m12 + latent SV"
)

mu_summary <- function(fit, nm, dates) {
  d <- unclass(as.matrix(fit$draws("mu", format = "draws_matrix")))
  tibble(
    date = dates, model = factor(labels_m[nm], levels = labels_m),
    mean = colMeans(d),
    lo = apply(d, 2, quantile, 0.05),
    hi = apply(d, 2, quantile, 0.95)
  )
}

mu_full <- bind_rows(lapply(names_full, function(nm) {
  mu_summary(readRDS(fit_path(nm)), nm, gdp_growth$date)
}))

p1 <- ggplot(mu_full, aes(x = date)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(data = gdp_growth, aes(y = headline),
            colour = "grey50", alpha = 0.6) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.25, fill = "steelblue") +
  geom_line(aes(y = mean), colour = "steelblue", linewidth = 0.4) +
  facet_wrap(~ model, ncol = 1) +
  labs(
    title = "Latent GDP growth (blue, 90% CI) vs headline (grey)",
    subtitle = paste0("Annualised quarterly growth, %. Headline = mean of GDP(E,I,P).", samp_note),
    x = NULL, y = NULL
  )
ggsave(file.path(fig_dir, "01_latent_vs_headline.png"),
       p1, width = 10, height = 17, dpi = 150)

# SV paths: m4 (independent errors), m9 (structural Sigma), m12 (m9 + AR(1) errors;
# its sigma_t is the sd of the error innovation u_t) and m13 (m12 + latent SV; same sigma_t), posterior mean and 90% interval of
# the measurement-error sd sigma_t (filename kept from the m4-only version of this figure).
T_ <- nrow(gdp_growth)
sv_df <- bind_rows(lapply(c("m4", "m9", "m12", "m13"), function(nm) {
  d <- unclass(as.matrix(readRDS(fit_path(nm))$draws("sigma_t", format = "draws_matrix")))
  stopifnot(ncol(d) == T_ * 3)                     # sigma_t[t, k], column-major
  tibble(
    model = factor(labels_m[nm], levels = labels_m),
    date = rep(gdp_growth$date, 3),
    approach = rep(c("E", "I", "P"), each = T_),
    mean = colMeans(d),
    lo = apply(d, 2, quantile, 0.05),
    hi = apply(d, 2, quantile, 0.95)
  )
}))
p2 <- ggplot(sv_df, aes(x = date, colour = approach, fill = approach)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(aes(y = mean), linewidth = 0.6) +
  facet_wrap(~ model, ncol = 1) +
  labs(title = "Measurement-error sd over time (posterior mean, 90% interval)",
       subtitle = paste0("Annualised %. SV models m4 (independent errors), m9 (structural Sigma), m12 (m9 + AR(1) errors: sd of the innovation) and m13 (m12 + latent SV).", samp_note),
       x = NULL, y = expression(sigma[t]), colour = "Approach", fill = "Approach")
ggsave(file.path(fig_dir, "02_m4_sv_path.png"),
       p2, width = 9, height = 12, dpi = 150)

# m7 kappa posterior
fit_m7 <- readRDS(fit_path("m7"))
kappa_draws <- fit_m7$draws("kappa", format = "draws_matrix")
covid_labels <- c("2020 Q2", "2020 Q3", "2020 Q4", "2021 Q1")
kappa_df <- as_tibble(unclass(as.matrix(kappa_draws))) |>
  setNames(covid_labels) |>
  pivot_longer(everything(), names_to = "quarter", values_to = "kappa") |>
  mutate(quarter = factor(quarter, levels = covid_labels))
p3 <- ggplot(kappa_df, aes(x = kappa, fill = quarter)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 1, colour = "grey50", linetype = 2) +
  facet_wrap(~ quarter, ncol = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 20)) +
  guides(fill = "none") +
  labs(title = "m7: posterior of COVID-quarter variance multipliers (kappa)",
       subtitle = paste0("Dashed line = no inflation. Only 2020 Q2 needed scale-up.", samp_note),
       x = expression(kappa), y = NULL)
ggsave(file.path(fig_dir, "04_m7_kappa.png"), p3, width = 8, height = 6, dpi = 150)

# 4-signal models (hours subsample): m8 and m10 latent paths vs the observed series
# (filename kept from the m8-only version of this figure).
mu_4 <- bind_rows(lapply(c("m8", "m10"), function(nm)
  mu_summary(readRDS(fit_path(nm)), nm, gdp_growth_4$date)))
gdp4_long <- gdp_growth_4 |>
  select(date, `GDP(E)` = g_E, `GDP(I)` = g_I, `GDP(P)` = g_P,
         `Hours worked` = g_U) |>
  pivot_longer(-date, names_to = "series", values_to = "growth")
p4 <- ggplot(mu_4, aes(x = date)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(data = gdp4_long, aes(y = growth, colour = series), alpha = 0.55) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "black") +
  geom_line(aes(y = mean), colour = "black", linewidth = 0.7) +
  facet_wrap(~ model, ncol = 1) +
  scale_colour_manual(values = c(
    "GDP(E)" = "#d95f02", "GDP(I)" = "#1b9e77",
    "GDP(P)" = "#7570b3", "Hours worked" = "#e7298a"
  )) +
  labs(
    title = "4-signal latent growth (black, 90% CI) vs observed series",
    subtitle = sprintf("Hours subsample from %s. Hours loaded as alpha + beta*mu (annualised %%).",
                       format(min(gdp_growth_4$date), "%Y")),
    x = NULL, y = "%", colour = NULL
  )
ggsave(file.path(fig_dir, "05_m8_4signal.png"), p4,
       width = 10, height = 8, dpi = 150)

# m5 COVID zoom (unchanged from before, but rebuilt for consistency)
mu_m5 <- mu_full |>
  filter(model == labels_m["m5"], date >= as.Date("2018-01-01"))
gdp_zoom <- gdp_growth |> filter(date >= as.Date("2018-01-01"))
p5 <- ggplot(mu_m5, aes(x = date)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(data = gdp_zoom, aes(y = g_E, colour = "GDP(E)"), alpha = 0.6) +
  geom_line(data = gdp_zoom, aes(y = g_I, colour = "GDP(I)"), alpha = 0.6) +
  geom_line(data = gdp_zoom, aes(y = g_P, colour = "GDP(P)"), alpha = 0.6) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "black") +
  geom_line(aes(y = mean), colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = c("GDP(E)" = "#d95f02",
                                  "GDP(I)" = "#1b9e77", "GDP(P)" = "#7570b3")) +
  labs(title = "m5 latent (black, 90% CI) and three GDP measures, 2018-onwards",
       x = NULL, y = "annualised %", colour = NULL)
ggsave(file.path(fig_dir, "03_m5_covid_zoom.png"), p5,
       width = 9, height = 5, dpi = 150)

# m11 monthly latent path
fit_m11 <- readRDS(fit_path("m11"))
emp_growth <- sw_trim_m(readRDS("data/processed/employment_growth_monthly.rds"))
gdp_growth_full <- sw_trim_q(readRDS("data/processed/gdp_growth.rds")) |>
  mutate(headline = (g_E + g_I + g_P) / 3)

mu_m_draws <- unclass(as.matrix(fit_m11$draws("mu_m", format = "draws_matrix")))
m_dates <- emp_growth |>
  filter(date <= max(gdp_growth_full$date)) |>
  pull(date)
stopifnot(length(m_dates) == ncol(mu_m_draws))   # fit_m11.rds must match current data
mu_m_df <- tibble(
  date = m_dates,
  mean = colMeans(mu_m_draws),
  lo   = apply(mu_m_draws, 2, quantile, 0.05),
  hi   = apply(mu_m_draws, 2, quantile, 0.95)
)
p6 <- ggplot(mu_m_df, aes(x = date)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(data = filter(gdp_growth_full, date >= min(m_dates)),
            aes(y = headline), colour = "grey55", alpha = 0.7) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.3, fill = "darkorange") +
  geom_line(aes(y = mean), colour = "darkorange", linewidth = 0.5) +
  labs(
    title = "m11: monthly latent growth (orange, 90% CI) vs quarterly headline (grey)",
    subtitle = paste0("Mixed-frequency model with monthly employment as a 4th signal", samp_note),
    x = NULL, y = "annualised %"
  )
ggsave(file.path(fig_dir, "06_m11_monthly_path.png"), p6,
       width = 10, height = 5, dpi = 150)

# PSIS-LFO: cumulative one-step-ahead log predictive density relative to a reference
# model, from R/03_compare.R (<root>/tables/lfo_pointwise.csv).
lfo_pw_path <- sw_out("tables", "lfo_pointwise.csv")
if (file.exists(lfo_pw_path)) {
  pw <- read.csv(lfo_pw_path) |> mutate(date = as.Date(date))
  rel <- function(d, ref, col = "elpd") {
    r <- d |> filter(model == ref) |> select(date, ref_v = all_of(col))
    d |> select(model, date, v = all_of(col)) |> inner_join(r, by = "date") |>
      group_by(model) |> arrange(date) |> mutate(cum = cumsum(v - ref_v)) |> ungroup()
  }
  g1 <- c("m1", "m2", "m3", "m4", "m5", "m6", "m7", "m9", "m12", "m13")
  d1 <- rel(filter(pw, model %in% g1), "m2") |> mutate(panel = "Group 1 (E, I, P), relative to m2")
  d2 <- rel(filter(pw, model %in% c("m8", "m10")), "m8") |>
    mutate(panel = "Group 2 (E, I, P, hours), relative to m8")
  # GDP-only target on the common window, relative to m12 (production): group-1 models'
  # joint density is the GDP density; m10's is its GDP-only column (hours in the past only).
  pg <- pw |> filter(model %in% c("m9", "m12", "m13", "m10")) |>
    mutate(elpd = ifelse(model == "m10", elpd_sub, elpd))
  d3 <- rel(pg, "m12") |> mutate(panel = "GDP-only target, relative to m12 (m10 also sees past hours)")
  panels <- c("Group 1 (E, I, P), relative to m2", "Group 2 (E, I, P, hours), relative to m8",
              "GDP-only target, relative to m12 (m10 also sees past hours)")
  dd <- bind_rows(d1, d2, d3) |>
    mutate(model = factor(labels_m[model], levels = labels_m), panel = factor(panel, levels = panels))
  # One fixed, clearly separated colour per model, so that the legend alone identifies
  # every line (ggplot's default hue scale gave m10 and m12 near-identical pinks).
  pal7 <- c(m1 = "#9e9e9e", m2 = "#000000", m3 = "#a6761d", m4 = "#56b4e9",
            m5 = "#e69f00", m6 = "#cc79a7", m7 = "#009e73", m8 = "#7f7f7f",
            m9 = "#0072b2", m10 = "#d55e00", m12 = "#7b3294", m13 = "#76c000")
  pal7 <- setNames(pal7, labels_m[names(pal7)])
  p7 <- ggplot(dd, aes(date, cum, colour = model)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_vline(xintercept = as.Date("2020-06-01"), colour = "grey80", linetype = 2) +
    geom_step(aes(linewidth = model == labels_m[["m13"]])) +
    scale_linewidth_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.2), guide = "none") +
    scale_colour_manual(values = pal7, drop = TRUE) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    labs(title = "PSIS-LFO: cumulative one-step-ahead log predictive density",
         subtitle = paste0("Marginal likelihood (latent states integrated out; SV models: mean of 4 particle-filter passes). Above 0 = better than the reference model.",
                           if (nzchar(samp_note)) paste0("\n", trimws(samp_note), " Targets 1990 Q1 onward, as on the full sample.")),
         x = NULL, y = "cumulative elpd difference", colour = NULL)
  ggsave(file.path(fig_dir, "07_lfo_cumulative_elpd.png"), p7, width = 10, height = 11, dpi = 150)
}

cat("Figures written to ", fig_dir, "/:\n",
    "  01_latent_vs_headline.png  (all full-sample specs including m9, m12 and m13)\n",
    "  02_m4_sv_path.png          (m4, m9, m12 and m13 measurement-error sd paths)\n",
    "  03_m5_covid_zoom.png\n",
    "  04_m7_kappa.png\n",
    "  05_m8_4signal.png          (m8 and m10)\n",
    "  06_m11_monthly_path.png\n",
    "  07_lfo_cumulative_elpd.png (if R/03_compare.R has run)\n", sep = "")
