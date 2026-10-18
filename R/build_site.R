# Generates the static site under docs/ (or $LATENTGDP_DOCS):
#   <docs>/index.html          - landing page with headline + chart
#   <docs>/chart_latent.png    - historical mu vs headline
#   <docs>/chart_recent.png    - zoomed view of last ~5 years
# Assumes output_estimate.R has already written <docs>/estimates*.json/csv,
# <docs>/diagnostics.json and <docs>/model_summary.json.
#
# Environment override (for smoke tests that must not touch docs/):
#   LATENTGDP_DOCS  site directory (default docs)

suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
  library(ggplot2)
})

theme_set(theme_minimal(base_size = 11))

docs_dir <- Sys.getenv("LATENTGDP_DOCS", "docs")

estimates <- read.csv(file.path(docs_dir, "estimates.csv")) |>
  mutate(date = as.Date(date))
latest    <- read_json(file.path(docs_dir, "estimates_latest.json"), simplifyVector = TRUE)
diag      <- read_json(file.path(docs_dir, "diagnostics.json"), simplifyVector = TRUE)
# Posterior summaries quoted in the model description (written by output_estimate.R
# from the production fit, so the text follows each new fit).
msum_path <- file.path(docs_dir, "model_summary.json")
if (!file.exists(msum_path))
  stop(msum_path, " not found: run R/output_estimate.R first", call. = FALSE)
msum      <- read_json(msum_path, simplifyVector = TRUE)

# "YYYY Qn" label for a quarter-end date.
quarter_label <- function(d) {
  d <- as.Date(d)
  sprintf("%s Q%d", format(d, "%Y"), (as.integer(format(d, "%m")) - 1) %/% 3 + 1)
}

# First quarter of the sample, for chart titles ("1959 Q4-present").
sample_start_q <- quarter_label(min(estimates$date))

# Y-axis range covering the 1960s-70s booms and normal recessions but
# clipping the COVID collapse; the annotation flags what's off-scale. Units are
# quarter-on-quarter growth (100 x dlog), 1/4 of the annualised scale.
YLIM <- c(-2.5, 4.5)

# Headline quarters outside YLIM are found from the data (not hard-coded),
# marked with a triangle at the plot edge and listed in an annotation (the
# most extreme few, if there are many).
off_scale <- estimates |>
  filter(headline_mean < YLIM[1] | headline_mean > YLIM[2]) |>
  mutate(y_edge = pmin(pmax(headline_mean, YLIM[1] + 0.08), YLIM[2] - 0.08),
         dir    = ifelse(headline_mean > YLIM[2], "up", "down"))
off_scale_note <- if (nrow(off_scale) > 0) {
  MAX_LISTED <- 4
  listed <- off_scale |>
    arrange(desc(abs(headline_mean))) |>
    head(MAX_LISTED) |>
    arrange(date)
  paste0(
    "Off-scale headline: ",
    paste(sprintf("%s %+.1f%%", quarter_label(listed$date),
                  listed$headline_mean), collapse = ", "),
    if (nrow(off_scale) > MAX_LISTED)
      sprintf(" (+%d more)", nrow(off_scale) - MAX_LISTED) else ""
  )
} else NULL

# Full sample chart
p_full <- ggplot(estimates, aes(x = date)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(aes(y = headline_mean), colour = "grey55", alpha = 0.7) +
  geom_ribbon(aes(ymin = mu_q05, ymax = mu_q95),
              alpha = 0.25, fill = "steelblue") +
  geom_line(aes(y = mu_median), colour = "steelblue", linewidth = 0.5) +
  coord_cartesian(ylim = YLIM) +
  {if (nrow(off_scale) > 0)
    geom_point(data = off_scale, aes(y = y_edge, shape = dir),
               colour = "grey40", size = 1.6, show.legend = FALSE)} +
  scale_shape_manual(values = c(up = 24, down = 25)) +
  {if (!is.null(off_scale_note))
    annotate("text", x = max(estimates$date), y = YLIM[2] - 0.25,
             label = off_scale_note, size = 3, hjust = 1, colour = "grey40")} +
  labs(
    title = sprintf("Australian latent GDP growth, %s-present", sample_start_q),
    subtitle = "Model m13 (structural Sigma + SV + AR(1) errors + latent SV). Blue = posterior median + 90% CI. Grey = headline.",
    x = NULL, y = "quarter-on-quarter growth, %",
    caption = paste("Data: ABS 5206.0 Table 24. Run:", diag$run_time_utc, "UTC.")
  )
ggsave(file.path(docs_dir, "chart_latent.png"), p_full, width = 10, height = 5, dpi = 150)

# Recent view: skip COVID entirely so the chart shows post-pandemic dynamics
# at a useful y-scale.
recent <- estimates |> filter(date >= as.Date("2022-01-01"))
p_recent <- ggplot(recent, aes(x = date)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_line(aes(y = g_E, colour = "GDP(E)"), alpha = 0.5) +
  geom_line(aes(y = g_I, colour = "GDP(I)"), alpha = 0.5) +
  geom_line(aes(y = g_P, colour = "GDP(P)"), alpha = 0.5) +
  geom_ribbon(aes(ymin = mu_q05, ymax = mu_q95),
              alpha = 0.3, fill = "black") +
  geom_line(aes(y = mu_median), colour = "black", linewidth = 0.8) +
  scale_colour_manual(values = c(
    "GDP(E)" = "#d95f02", "GDP(I)" = "#1b9e77", "GDP(P)" = "#7570b3"
  )) +
  labs(
    title = "Latent GDP vs three ABS measures - post-2022",
    x = NULL, y = "quarter-on-quarter %", colour = NULL
  )
ggsave(file.path(docs_dir, "chart_recent.png"), p_recent, width = 10, height = 5, dpi = 150)

# Pretty-print the headline number with sign-aware colour for HTML; an
# estimate that rounds to 0.00 is shown unsigned in a neutral colour.
mu <- latest$mu_median
mu_colour <- if (abs(mu) < 0.005) "#555555" else if (mu > 0) "#1b9e77" else "#d62728"
mu_str <- if (abs(mu) < 0.005) "0.00%" else sprintf("%+.2f%%", mu)
ci_str <- sprintf("[%+.2f, %+.2f]", latest$mu_q05, latest$mu_q95)

# Variance of the latent path relative to the headline: posterior median of
# the per-draw ratio, with its 90% interval when available.
pct <- function(x) sprintf("%.0f%%", 100 * x)
var_str <- if (!is.null(latest$variance_ratio_vs_headline_q05) &&
               !is.null(latest$variance_ratio_vs_headline_q95)) {
  sprintf("%s of headline variance across the full sample (posterior median; 90%% interval %s&ndash;%s)",
          pct(latest$variance_ratio_vs_headline),
          pct(latest$variance_ratio_vs_headline_q05),
          pct(latest$variance_ratio_vs_headline_q95))
} else {
  sprintf("%s of headline variance across the full sample",
          pct(latest$variance_ratio_vs_headline))
}

# Footer diagnostics, read defensively (older diagnostics.json lacked fields).
fmt_or_na <- function(x, fmt) if (is.null(x) || length(x) == 0 || is.na(x)) "n/a" else sprintf(fmt, x)
n_chains    <- if (!is.null(diag$chains)) diag$chains else NA
n_per_chain <- if (!is.null(diag$draws_per_chain)) diag$draws_per_chain else diag$iter_per_chain
total_draws <- if (!is.null(diag$total_draws)) diag$total_draws else n_chains * n_per_chain
n_div <- as.integer(diag$num_divergent)
footer_diag <- sprintf(
  "%d %s out of %s post-warmup draws (%s chains &times; %s). Min E-BFMI %s; max R-hat (hyperparameters) %s; min bulk / tail ESS %s / %s.",
  n_div, if (identical(n_div, 1L)) "divergence" else "divergences",
  fmt_or_na(total_draws, "%d"),
  fmt_or_na(n_chains, "%d"), fmt_or_na(n_per_chain, "%d"),
  fmt_or_na(diag$ebfmi_min, "%.3f"), fmt_or_na(diag$rhat_max_hyper, "%.3f"),
  fmt_or_na(diag$ess_bulk_min, "%.0f"), fmt_or_na(diag$ess_tail_min, "%.0f")
)

# The page is assembled from three parts, because R's sprintf() limits a format
# string to 8192 bytes: html_top and html_foot are sprintf() templates (so a
# literal percent sign is written %%), html_model is plain text (a single %).
css <- '
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           max-width: 880px; margin: 2rem auto; padding: 0 1rem; color: #222;
           line-height: 1.5; }
    h1 { font-size: 1.5rem; margin-bottom: 0.5rem; }
    .headline { font-size: 3rem; font-weight: 600; margin: 0; }
    .ci { color: #666; font-size: 1.1rem; }
    .meta { color: #888; font-size: 0.9rem; margin-top: 0.5rem; }
    .card { background: #f7f7f7; padding: 1.2rem 1.5rem;
            border-radius: 8px; margin: 1rem 0; }
    table { border-collapse: collapse; margin: 0.5rem 0; }
    td { padding: 0.2rem 0.6rem 0.2rem 0; }
    td:first-child { color: #666; }
    img { max-width: 100%; height: auto; display: block; margin: 1rem 0; }
    code { background: #eee; padding: 0.1rem 0.3rem; border-radius: 3px; }
    .eq { background: #f7f7f7; border-left: 3px solid #c0c8d8;
          padding: 0.5rem 0.9rem; margin: 0.5rem 0;
          font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
          font-size: 0.92rem; white-space: pre-wrap; overflow-x: auto; }
    a { color: #1f6feb; }
    footer { color: #888; font-size: 0.85rem; margin-top: 2rem;
             border-top: 1px solid #eee; padding-top: 1rem; }
'

html_top <- sprintf('<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Latent Australian GDP growth</title>
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <style>%s</style>
</head>
<body>
  <h1>Latent Australian GDP growth, %s</h1>
  <div class="card">
    <p class="headline" style="color: %s;">%s</p>
    <p class="ci">90%% credible interval: %s (quarter-on-quarter growth)</p>
    <p class="meta">Headline (mean of E, I, P): %+.2f%%. Latent-path variance is %s.</p>
  </div>

  <h2>Observed measures for %s</h2>
  <table>
    <tr><td>GDP(E) expenditure-based</td><td><strong>%+.2f%%</strong></td></tr>
    <tr><td>GDP(I) income-based</td>     <td><strong>%+.2f%%</strong></td></tr>
    <tr><td>GDP(P) production-based</td> <td><strong>%+.2f%%</strong></td></tr>
  </table>

  <h2>Recent quarters (post-2022)</h2>
  <img src="chart_recent.png" alt="Latent GDP growth vs three measures, post-2022">

  <h2>About</h2>
  <p>This page publishes a Bayesian state-space estimate of Australian
  quarterly GDP growth that filters measurement error out of the three
  ABS measures (Expenditure, Income, Production-based GDP).
  The model is m13 from <a href="https://github.com/DavidAStephan/LatentGDP">this
  repo</a>: an AR(1) latent state whose volatility changes over time,
  measurement errors that follow their own AR(1) with stochastic volatility,
  and a construction-based zero restriction on the Production-side error
  covariance. Background: RBA RDP 2014-12 (Rees, Lancaster &amp; Finlay).</p>

  <p>Data: ABS Cat. 5206.0 Table 24 (Selected Analytical Series),
  seasonally adjusted, chain volume measures. The page is rebuilt
  automatically after each quarterly release.</p>

  <h3>Historical estimate, %s&ndash;present</h3>
  <img src="chart_latent.png" alt="Latent GDP growth vs headline, %s to present">

  <h3>Downloads</h3>
  <ul>
    <li><a href="estimates.csv">estimates.csv</a> - full quarterly time series with 90%% CI</li>
    <li><a href="estimates_latest.json">estimates_latest.json</a> - latest quarter headline</li>
    <li><a href="diagnostics.json">diagnostics.json</a> - sampler diagnostics for the latest fit</li>
  </ul>
',
  css, latest$quarter, mu_colour, mu_str, ci_str,
  latest$headline_mean, var_str,
  latest$quarter,
  latest$observed$g_E, latest$observed$g_I, latest$observed$g_P,
  sample_start_q, sample_start_q
)

html_model <- '
  <h2>Model details</h2>

  <h3>Observation equations</h3>
  <div class="eq">y<sub>t</sub><sup>k</sup> = &mu;<sub>t</sub> + &epsilon;<sub>t</sub><sup>k</sup>,  k = E, I, P
&epsilon;<sub>t</sub><sup>k</sup> = &psi;<sub>k</sub> &epsilon;<sub>t-1</sub><sup>k</sup> + u<sub>t</sub><sup>k</sup>
(u<sub>t</sub><sup>E</sup>, u<sub>t</sub><sup>I</sup>)  ~  MVN(0, &Sigma;<sub>EI,t</sub>),   u<sub>t</sub><sup>P</sup>  ~  N(0, &sigma;<sub>P,t</sub><sup>2</sup>)</div>
  <p>Each measurement error follows an AR(1) with coefficient
  &psi;<sub>k</sub>, so an error in one quarter can be partly carried into,
  or reversed in, the next. &Sigma;<sub>EI,t</sub> has standard deviations
  &sigma;<sub>E,t</sub>, &sigma;<sub>I,t</sub> and a constant correlation
  &rho;<sub>EI</sub>: E and I share national-accounts balancing items. P is
  built from independent industry surveys, so its error is independent of
  E and I by construction (the &quot;structural Sigma&quot;). The first
  quarter&apos;s errors start from the stationary distribution at that
  quarter&apos;s variances.</p>

  <h3>Latent state and stochastic volatility</h3>
  <div class="eq">&mu;<sub>t</sub> = c + &rho; &middot; &mu;<sub>t-1</sub> + &tau;<sub>t</sub> &eta;<sub>t</sub>,  &eta;<sub>t</sub> ~ N(0, 1),  &tau;<sub>t</sub> = exp(g<sub>t</sub> / 2)
&mu;<sub>1</sub> ~ N(c / (1 - &rho;), &tau;<sub>1</sub><sup>2</sup> / (1 - &rho;<sup>2</sup>))
g<sub>t</sub> = &mu;<sub>g</sub> + &phi;<sub>g</sub> (g<sub>t-1</sub> - &mu;<sub>g</sub>) + &omega;<sub>g</sub> &zeta;<sub>t</sub>,  &zeta;<sub>t</sub> ~ N(0, 1),  g<sub>1</sub> stationary
h<sub>k,t</sub> = &mu;<sub>h,k</sub> + &phi;<sub>k</sub> (h<sub>k,t-1</sub> - &mu;<sub>h,k</sub>) + &omega;<sub>k</sub> &xi;<sub>k,t</sub>,  &xi;<sub>k,t</sub> ~ N(0, 1),  h<sub>k,1</sub> stationary
&sigma;<sub>k,t</sub><sup>2</sup> = exp(h<sub>k,t</sub>) = Var(u<sub>t</sub><sup>k</sup>),  k = E, I, P</div>
  <p>There are two kinds of volatility. The measurement-error volatility
  h captures the long fall in measurement error (largest in the
  1960s&ndash;70s). The latent volatility &tau;<sub>t</sub> lets the
  volatility of true growth change too: in the current fit it
  falls from about {{TAU_1960S}} (annualised, decade average) in the 1960s to about
  {{TAU_2010S}} in the 2010s and rises to about {{TAU_2020Q2}} in 2020 Q2. COVID is therefore treated as a
  large real shock, not as measurement error: the three measures agreed in
  2020 Q2, so the collapse goes into &mu;<sub>t</sub>.
  The estimated &psi; are about {{PSI_E}}, {{PSI_I}} and {{PSI_P}} for E, I and P
  (posterior medians); negative values mean that the
  measurement errors partly reverse the next quarter.
  Given both volatility paths, &mu; is integrated out by a Kalman filter on
  the quasi-differenced data y<sub>t</sub> &minus; &psi; y<sub>t-1</sub>
  during sampling, and drawn by forward-filtering backward-sampling.
  With &tau;<sub>t</sub> held constant the model is m12, the previous
  production model.</p>

  <h3>Priors</h3>
  <p>In the model&apos;s annualised units (400&times;dlog; this page shows
  quarter-on-quarter, i.e. divided by 4). N(m, s) has sd s.</p>
  <div class="eq">c           ~ N(1.5, 1)
&rho;           ~ N(0.5, 0.3)    constrained to (-1, 1)
&mu;<sub>g</sub>          ~ N(2, 1)
&phi;<sub>g</sub>          ~ Uniform(-1, 1)
&omega;<sub>g</sub>          ~ half-N(0, 0.5)
&mu;<sub>h,k</sub>         ~ N(2, 1)
&phi;<sub>k</sub>          ~ Uniform(-1, 1)
&omega;<sub>k</sub>          ~ half-N(0, 0.5)
&rho;<sub>EI</sub>         ~ N(0, 0.5)    constrained to (-1, 1)
&psi;<sub>k</sub>          ~ N(0, 0.5)    constrained to (-1, 1)</div>
  <p>All priors for all models:
  <a href="https://github.com/DavidAStephan/LatentGDP/blob/master/priors.md">priors.md</a>.
  Nothing is published unless the fit has &lt;0.1% divergences, E-BFMI
  &ge; 0.3, R-hat &le; 1.01 (hyperparameters) / 1.05 (latent growth and
  both volatility paths) and ESS &ge; 400.</p>

  <h3>Why this specification</h3>
  <p>The 13 specifications were compared by leave-future-out
  cross-validation: one-step-ahead predictive density, latent states
  integrated out, 146 quarters from 1990 Q1 (current-vintage data). They
  were compared in three groups with different data: the three measures
  on the full sample (m1&ndash;m7, m9, m12, m13); the three measures plus
  hours from 1978 (m8, m10); and the mixed-frequency m11, which also sees
  monthly employment inside the predicted quarter.</p>
  <ul>
    <li>Three features help, each measured as one step added to a model
      that has the ones before (whether the gains add up was not tested):
      stochastic volatility
      in the measurement errors, serially correlated measurement errors,
      and stochastic volatility in latent growth. m13 has all three. It
      beats m12 (the first two only) by 59.6 nats (SE 24.7 assuming
      independent quarters, 28.9 allowing for serial correlation; z = 2.1).
      One quarter, 2020 Q2, contributes 23.7 nats and most of that SE.
      Without the COVID quarters the gain is 28.5 nats (SE 6.0), and m13 is
      ahead in 122 of the 146 quarters.</li>
    <li>m13 replaced m12 in October 2026 under a rule fixed before
      m13&apos;s comparison results were known: pass the publication checks
      on two seeds, have a reliable score, beat m12 by at least 2
      serial-correlation-robust SEs (also when no refit is allowed a retry
      with a new seed), and fit in under 30 minutes. It met every
      condition, with thin margins (z = 2.06; 23 minutes). m13 was proposed after the hours
      model m10, re-specified with the same latent volatility, had been
      scored on the same quarters, so its gain is not a pre-registered
      test. m12 had replaced m9 the day before under the same kind of rule
      (+16.8 nats, z = 2.7).</li>
    <li>The structural restriction itself is not supported by prediction
      (m9 against m4, the same model without it: +1.5 nats, SE 1.5 assuming
      independent quarters, 2.2 allowing for serial correlation). It is
      kept because it reflects how GDP(P) is built.</li>
    <li>On GDP prediction, the hours model m10, which also has latent
      volatility, is behind m13 by 12.6 nats (SE 5.5) and the employment
      model m11 is behind by 68 (SE 18; both SEs allow for serial
      correlation). These gaps do not show that hours or employment are
      unhelpful: m10 lacks m13&apos;s serially correlated measurement errors
      and E/I-only error correlation and is estimated only from 1978, m11
      has no stochastic volatility, and no m13-plus-hours model has been
      fitted.</li>
  </ul>
  <p>The comparison scores predictions of the observed measures, not the
  split between true growth and measurement error, and models that
  predict about equally well can imply very different splits. Over the full sample
  m13&apos;s latent path has {{VR_FULL}}% of headline variance (from the
  current production fit) and m12&apos;s 78%, but over 1978&ndash;2026
  excluding 2020 Q1&ndash;2021 Q1 m13&apos;s has 75%, m12&apos;s 88% and the
  hours model m10&apos;s 56%. Apart from m13&apos;s full-sample figure, these
  splits come from research-suite fits on the 2026 Q2 data and are not
  updated by the weekly run. The variance figure above is
  therefore a property of m13, not something the comparison confirms.</p>
  <p>An earlier PSIS-LOO comparison on this page (m9 +6 over m4, SV
  +240) was invalid and has been withdrawn.</p>

  <p>Full research record (all 13 specs, the comparison tables and the
  caveats) at
  <a href="https://github.com/DavidAStephan/LatentGDP/blob/master/NOTES.md">NOTES.md</a>.
  Stan source at
  <a href="https://github.com/DavidAStephan/LatentGDP/blob/master/stan/m13_serial_latent_sv.stan">stan/m13_serial_latent_sv.stan</a>.</p>
'

html_foot <- sprintf('
  <footer>
    Last updated %s UTC. CmdStan %s. %s
    Sample %s to %s (T = %d quarters). Source: <a href="https://github.com/DavidAStephan/LatentGDP">github.com/DavidAStephan/LatentGDP</a>.
  </footer>
</body>
</html>
',
  diag$run_time_utc, diag$cmdstan_version, footer_diag,
  quarter_label(diag$sample_start), quarter_label(diag$sample_end), as.integer(diag$T)
)

# Repository links point at the commit being built when run by GitHub Actions
# (GITHUB_SHA), so they resolve whichever branch the workflow runs on (e.g. a
# workflow_dispatch on a branch whose files are not yet on master); locally they
# point at master.
repo_ref <- Sys.getenv("GITHUB_SHA", "")
if (!grepl("^[0-9a-f]{40}$", repo_ref)) repo_ref <- "master"
html_model <- gsub("/LatentGDP/blob/master/", sprintf("/LatentGDP/blob/%s/", repo_ref),
                   html_model, fixed = TRUE)

# Fill the posterior numbers in the model description from model_summary.json.
# (The variance splits of other models in "Why this specification" are research-suite
# results, labelled as such, and are edited by hand when the research suite is rerun.)
num <- function(x, fmt) {
  if (is.null(x) || length(x) != 1 || !is.finite(x)) stop("model_summary.json: missing value", call. = FALSE)
  s <- sprintf(fmt, x)
  if (x < 0) paste0("&minus;", sub("^-", "", s)) else s
}
fills <- c(
  "{{TAU_1960S}}"  = num(msum$latent_sd_decade_avg[["1960s"]], "%.1f"),
  "{{TAU_2010S}}"  = num(msum$latent_sd_decade_avg[["2010s"]], "%.1f"),
  "{{TAU_2020Q2}}" = num(msum$latent_sd_2020q2, "%.0f"),
  "{{PSI_E}}"      = num(msum$psi[["E"]], "%.2f"),
  "{{PSI_I}}"      = num(msum$psi[["I"]], "%.2f"),
  "{{PSI_P}}"      = num(msum$psi[["P"]], "%.2f"),
  "{{VR_FULL}}"    = num(100 * latest$variance_ratio_vs_headline, "%.0f"))
for (k in names(fills)) html_model <- gsub(k, fills[[k]], html_model, fixed = TRUE)
stopifnot(!grepl("{{", html_model, fixed = TRUE))

html <- paste0(html_top, html_model, html_foot)

writeLines(html, file.path(docs_dir, "index.html"))
cat("Site built at ", docs_dir, "/:\n",
    "  index.html, chart_latent.png, chart_recent.png\n",
    "  estimates.csv, estimates_latest.json, diagnostics.json, model_summary.json\n", sep = "")
