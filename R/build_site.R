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

# First quarter of the sample, for chart titles ("1980 Q1-present").
sample_start_q <- quarter_label(min(estimates$date))

# Y-axis range covering normal recessions and the 2020 Q3 rebound but
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
    subtitle = "Blue = estimated true growth (median and 90% interval). Grey = headline (average of the three ABS measures).",
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
  sprintf("%s of headline variance over %s&ndash;present (posterior median; 90%% interval %s&ndash;%s)",
          pct(latest$variance_ratio_vs_headline), sample_start_q,
          pct(latest$variance_ratio_vs_headline_q05),
          pct(latest$variance_ratio_vs_headline_q95))
} else {
  sprintf("%s of headline variance over %s&ndash;present",
          pct(latest$variance_ratio_vs_headline), sample_start_q)
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
  <p>This page publishes an estimate of Australian quarterly GDP growth
  with measurement error filtered out of the three ABS measures
  (expenditure, income and production based GDP). How the estimate is
  made is described below. Background: RBA Research Discussion Paper
  2014-12 (Rees, Lancaster &amp; Finlay).</p>

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
  <h2>How the estimate is made</h2>

  <p>Each quarter the ABS publishes three estimates of GDP growth built
  from different sources: spending (E), incomes (I) and industry
  production (P). They never agree exactly, because each contains
  measurement error. The model treats true GDP growth as an unobserved
  series, treats each published measure as true growth plus its own error,
  and works out the most likely path of true growth given all three.</p>

  <p>It rests on four assumptions:</p>
  <ul>
    <li><strong>True growth</strong> moves around a long-run average and
      partly persists from one quarter to the next. How volatile it is can
      change over time. The typical size of a quarterly shock to true growth
      is estimated at about {{TAU_1980S}} percentage points in the 1980s,
      {{TAU_2010S}} in the 2010s and {{TAU_2020Q2}} in 2020 Q2. Because the
      three measures agreed on the size of the 2020 Q2 fall, the model
      treats COVID as a genuine collapse in activity rather than a
      measurement problem.</li>
    <li><strong>Measurement errors change in size over time.</strong> Each
      measure&apos;s errors can be larger in some periods than in others.</li>
    {{PSI_ITEM}}
    <li><strong>Only some errors are linked.</strong> The spending and
      income measures share some of the same balancing adjustments, so their
      errors may move together. The production measure is built from
      separate industry surveys, so its errors are assumed independent of
      the other two.</li>
  </ul>

  <p>The model is estimated on data from 1980 Q1 onward, following the
  RBA paper, because measurement errors before then may not reflect how the
  ABS measures behave today. Over that period the estimated true-growth
  path has {{VR_FULL}}% of the variance of headline growth, so roughly
  {{NOISE}}% of the variance of quarterly headline growth is measurement
  noise. Estimation is Bayesian, so every figure comes with an
  interval that shows its uncertainty.</p>

  <p>Before anything is published, the estimate has to pass automatic
  checks that the statistical sampling worked properly. If any check
  fails, this page is not updated.</p>

  <p>The model was chosen by comparing a range of alternatives on how well
  each forecast the next quarter&apos;s published figures. The comparison,
  its caveats and the code are on
  <a href="https://github.com/DavidAStephan/LatentGDP/blob/master/NOTES.md">GitHub</a>.</p>
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

# Fill the posterior numbers in the model description from model_summary.json and
# the latest estimate, so the text always matches the published fit. Latent
# volatility is converted from the model's annualised units to quarter-on-quarter
# percentage points (divide by 4), the units used everywhere else on the page.
num <- function(x, fmt) {
  if (is.null(x) || length(x) != 1 || !is.finite(x)) stop("model_summary.json: missing value", call. = FALSE)
  s <- sprintf(fmt, x)
  if (x < 0) paste0("&minus;", sub("^-", "", s)) else s
}
# The list item on how measurement errors carry into the next quarter, written
# from the signs of the estimated AR(1) coefficients psi (E, I, P).
psi_item <- function(psi) {
  v <- vapply(c("E", "I", "P"), function(k) {
    x <- psi[[k]]
    if (is.null(x) || length(x) != 1 || !is.finite(x)) stop("model_summary.json: missing psi", call. = FALSE)
    x
  }, numeric(1))
  pc <- sprintf("%.0f%%", 100 * abs(v))
  if (all(v < 0)) {
    sprintf(paste0('<li><strong>Errors partly unwind.</strong> In the quarter after an error,\n',
                   '      about %s of it is reversed for the spending measure, %s for\n',
                   '      income and %s for production.</li>'), pc[1], pc[2], pc[3])
  } else {
    dir <- ifelse(v < 0, "reversed", "carried over")
    sprintf(paste0('<li><strong>Errors partly carry over.</strong> In the quarter after an\n',
                   '      error, about %s of it is %s for the spending measure, %s is %s\n',
                   '      for income and %s is %s for production.</li>'),
            pc[1], dir[1], pc[2], dir[2], pc[3], dir[3])
  }
}
vr <- latest$variance_ratio_vs_headline
fills <- c(
  "{{TAU_1980S}}"  = num(msum$latent_sd_decade_avg[["1980s"]] / 4, "%.1f"),
  "{{TAU_2010S}}"  = num(msum$latent_sd_decade_avg[["2010s"]] / 4, "%.1f"),
  "{{TAU_2020Q2}}" = num(msum$latent_sd_2020q2 / 4, "%.1f"),
  "{{PSI_ITEM}}"   = psi_item(msum$psi),
  "{{VR_FULL}}"    = num(100 * vr, "%.0f"),
  "{{NOISE}}"      = num(100 - round(100 * vr), "%.0f"))
for (k in names(fills)) html_model <- gsub(k, fills[[k]], html_model, fixed = TRUE)
stopifnot(!grepl("{{", html_model, fixed = TRUE))

html <- paste0(html_top, html_model, html_foot)

writeLines(html, file.path(docs_dir, "index.html"))
cat("Site built at ", docs_dir, "/:\n",
    "  index.html, chart_latent.png, chart_recent.png\n",
    "  estimates.csv, estimates_latest.json, diagnostics.json, model_summary.json\n", sep = "")
