# Shared helpers for the ABS fetch scripts (R/01_fetch_data.R and
# R/01b_fetch_research_data.R). Source from the repo root:
#   source("R/utils_abs.R")

# --- CI resilience -----------------------------------------------------------
# read_abs() runs check_abs_connection() before downloading anything. That
# pre-flight check pings the www.abs.gov.au *homepage* (httr HEAD/GET plus a
# base-R url() fallback) and aborts if it doesn't get a 2xx. From cloud/CI IPs
# the ABS WAF (CloudFront/Varnish) intermittently answers the homepage with a
# non-2xx code, so the check fails even though the actual data hosts
# (ausstats.abs.gov.au time-series directory + the file downloads) are
# reachable. This is exactly what broke the scheduled run on 2026-06-03.
# Neutralise the fragile homepage gate so a false negative there can't kill the
# whole pipeline, and wrap every fetch in retry-with-backoff so a genuinely
# transient WAF challenge on the data host self-heals instead of failing CI.
neutralise_abs_connection_check <- function() {
  suppressWarnings(
    assignInNamespace("check_abs_connection", function() invisible(TRUE),
                      "readabs")
  )
  invisible(TRUE)
}

with_retry <- function(fn, what, tries = 5, base_sleep = 15) {
  for (attempt in seq_len(tries)) {
    res <- tryCatch(fn(), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    msg <- conditionMessage(res)
    if (attempt == tries) {
      stop(sprintf("%s failed after %d attempts. Last error: %s",
                   what, tries, msg), call. = FALSE)
    }
    sleep <- base_sleep * attempt   # 15s, 30s, 45s, 60s
    message(sprintf("[retry] %s failed (attempt %d/%d): %s -- retrying in %ds",
                    what, attempt, tries, msg, sleep))
    Sys.sleep(sleep)
  }
}

# Quarterly log-difference, annualised (x400) - matches paper units.
annualised_qoq <- function(x) 400 * (log(x) - log(dplyr::lag(x)))
