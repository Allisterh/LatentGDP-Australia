# Production sample: the single place that defines which quarters the published
# estimate is fitted on. Sourced by R/fit_production.R, R/output_estimate.R and
# R/input_hash.R, so the fit, the published series and the input fingerprint
# always use the same data.
#
# The first growth observation is 1980 Q1 (date 1980-03-01, i.e.
# 400 x (log level 1980 Q1 - log level 1979 Q4)); the end is the latest release.
# This follows RBA RDP 2014-12 (Rees, Lancaster & Finlay), section 3.2: "Our data
# span 1980:Q1-2013:Q2. The starting date reflects the fact that, while Australian
# national accounts data are available on a quarterly basis from 1959:Q3, the
# quality of the underlying data sources has changed over time, so that the
# pattern of measurement errors in the early years of each GDP series may be
# unrepresentative of their current performance."
#
# Adopted 2026-10-19. The research suite's 1980 Q1 rerun (NOTES.md,
# output/sample_1980/) found that m13 passes the production gate on two seeds on
# this sample and still ranks first in the LFO comparison. The research suite
# itself (R/02_fit_models.R etc.) still defaults to the full sample; see
# R/sample_window.R for LATENTGDP_START.

PRODUCTION_START <- as.Date("1980-03-01")

load_production_growth <- function(path = "data/processed/gdp_growth.rds") {
  g <- readRDS(path)
  g <- g[as.Date(g$date) >= PRODUCTION_START, , drop = FALSE]
  if (!nrow(g) || as.Date(g$date[1]) != PRODUCTION_START) {
    stop("production_config: ", path, " has no observation dated ",
         format(PRODUCTION_START), call. = FALSE)
  }
  rownames(g) <- NULL
  g
}
