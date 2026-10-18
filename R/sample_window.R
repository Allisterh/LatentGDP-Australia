# Sample window and output root for the research suite (NOT used by production: the
# production scripts R/fit_production.R, R/output_estimate.R, R/build_site.R and
# R/input_hash.R never source this file and always use the full sample).
#
#   source("R/sample_window.R")
#
# Two environment variables, both optional:
#
#   LATENTGDP_START  ISO date of the FIRST GROWTH OBSERVATION, in the quarterly date
#                    convention of data/processed/gdp_growth.rds (first day of the last
#                    month of the quarter). "1980-03-01" makes 1980 Q1 the first
#                    observation, i.e. 400 x (log level 1980 Q1 - log level 1979 Q4); this
#                    is the sample start of RBA RDP 2014-12 (Rees, Lancaster & Finlay,
#                    section 3.2: "Our data span 1980:Q1-2013:Q2"). Unset or empty = the
#                    full sample (unchanged behaviour). The date must be one of the
#                    quarterly dates in gdp_growth.rds.
#   LATENTGDP_OUT    Output root (default "output"). Fits go to <root>/fit_*.rds, tables to
#                    <root>/tables/, figures to <root>/figures/, LFO checkpoints to
#                    <root>/lfo/. Use a separate root for any non-default sample, e.g.
#                    LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980.
#                    Sourcing this file stops with an error if LATENTGDP_START is set and
#                    the root is output/ (sw_guard() below).
#
# With both unset every function below returns its input unchanged (identical(), not
# just all.equal()) and sw_out("x") is "output/x", so the scripts that use it behave
# exactly as before.
#
# Truncation is the same for every data set: keep the quarters dated >= start.
#   - gdp_growth (m1-m7, m9, m12, m13) and gdp_growth_with_hours (m8, m10): sw_trim_q().
#     The hours sample starts in 1978 Q4; a start earlier than a series' first date
#     leaves that series unchanged (so the default full-sample hours subsample is kept).
#   - monthly employment (m11): sw_trim_m() keeps the months from the first month of the
#     start quarter (1980-01 for start 1980-03-01), so the start quarter is complete;
#     R/fit_m11.R then keeps only quarters whose three months are all present.
# The end of the sample is never changed (the latest data).
#
# Difference from the RDP: its fourth signal is the unemployment rate; m8 and m10 here
# use hours worked (5206.0 Table 1). Only the sample start is aligned with the paper.

sw_start <- function() {
  s <- Sys.getenv("LATENTGDP_START", "")
  if (!nzchar(s)) return(NULL)
  d <- tryCatch(as.Date(s, format = "%Y-%m-%d"), error = function(e) NA)
  if (is.na(d) || format(d, "%Y-%m-%d") != s)
    stop("LATENTGDP_START must be an ISO date (YYYY-MM-DD), got '", s, "'", call. = FALSE)
  if (!as.integer(format(d, "%m")) %in% c(3L, 6L, 9L, 12L) || format(d, "%d") != "01")
    stop("LATENTGDP_START = ", s, " is not a quarterly date: use the first day of the ",
         "quarter's last month (1980 Q1 = 1980-03-01)", call. = FALSE)
  d
}

sw_out_root <- function() {
  o <- Sys.getenv("LATENTGDP_OUT", "")
  if (!nzchar(o)) return("output")
  sub("/+$", "", o)
}

sw_is_default <- function() is.null(sw_start()) && sw_out_root() == "output"

# Path under the output root: sw_out("tables", "x.csv") = "<root>/tables/x.csv".
sw_out <- function(...) file.path(sw_out_root(), ...)

# Quarterly data frame with a `date` column: rows dated >= start.
sw_trim_q <- function(df, start = sw_start()) {
  if (is.null(start)) return(df)
  stopifnot(inherits(df$date, "Date"), !is.unsorted(df$date))
  if (start <= min(df$date)) return(df)
  if (!start %in% df$date)
    stop("LATENTGDP_START = ", format(start), " is not a date of this series (",
         format(min(df$date)), " to ", format(max(df$date)), ")", call. = FALSE)
  out <- df[df$date >= start, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Monthly data frame (date = first day of the month): months from the first month of the
# start quarter. The quarter's three months are then all inside the sample.
sw_trim_m <- function(df, start = sw_start()) {
  if (is.null(start)) return(df)
  stopifnot(inherits(df$date, "Date"), !is.unsorted(df$date))
  m0 <- seq(start, length.out = 2, by = "-2 months")[2]
  if (m0 <= min(df$date)) return(df)
  out <- df[df$date >= m0, , drop = FALSE]
  rownames(out) <- NULL
  out
}

# The three data sets of the research suite, trimmed consistently. Each element is read
# only if its file exists (hours and employment come from R/01b_fetch_research_data.R).
sw_data <- function(dir = "data/processed") {
  rd <- function(f) { p <- file.path(dir, f); if (file.exists(p)) readRDS(p) else NULL }
  g  <- rd("gdp_growth.rds")
  g4 <- rd("gdp_growth_with_hours.rds")
  em <- rd("employment_growth_monthly.rds")
  list(gdp_growth = if (!is.null(g)) sw_trim_q(g),
       gdp_growth_with_hours = if (!is.null(g4)) sw_trim_q(g4),
       employment_growth_monthly = if (!is.null(em)) sw_trim_m(em))
}

# Guard, run when this file is sourced (so it covers every research script that uses it:
# R/02_fit_models.R, R/fit_m11.R, R/fit_m12.R, R/fit_m13.R, R/03_compare.R, R/04_figures.R).
# A truncated sample must never write into the full-sample root output/: without this, a
# run with LATENTGDP_START set and LATENTGDP_OUT unset would overwrite output/fit_<m>.rds
# and its row of output/tables/fit_diagnostics.csv. The comparison is on normalized paths,
# so "output/", "./output" or the absolute path are caught too. It stops before any script
# has created a directory or written a file.
sw_guard <- function() {
  if (is.null(sw_start())) return(invisible(TRUE))
  root <- normalizePath(sw_out_root(), mustWork = FALSE)
  if (identical(root, normalizePath("output", mustWork = FALSE)))
    stop("LATENTGDP_START is set but the output root is output/ (LATENTGDP_OUT unset or ",
         "pointing at output/): a non-default sample must not write to output/ (use e.g. ",
         "LATENTGDP_OUT=output/sample_1980)", call. = FALSE)
  invisible(TRUE)
}
sw_guard()

# One line for logs and figure subtitles.
sw_describe <- function() {
  s <- sw_start()
  sprintf("sample start %s; output root %s",
          if (is.null(s)) "full sample (LATENTGDP_START unset)" else format(s), sw_out_root())
}
