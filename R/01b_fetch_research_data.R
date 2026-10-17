# Research-suite data fetch (NOT needed by the production pipeline / CI).
# Run AFTER R/01_fetch_data.R, which writes data/processed/gdp_growth.rds:
#   Rscript R/01_fetch_data.R
#   Rscript R/01b_fetch_research_data.R
#
# Writes:
#   data/processed/gdp_growth_with_hours.rds    - gdp_growth + g_U (hours
#       worked, 5206.0 Table 1), inner-joined on date; used by m8 / m10.
#   data/processed/employment_monthly.rds       - 6202.0 Table 1 employed
#       persons, SA level; used by m11.
#   data/processed/employment_growth_monthly.rds - month-on-month annualised
#       growth (1200 x dlog) of the above; used by m11.

suppressPackageStartupMessages({
  library(readabs)
  library(dplyr)
  library(tidyr)
})

source("R/utils_abs.R")
neutralise_abs_connection_check()

gdp_growth_path <- "data/processed/gdp_growth.rds"
if (!file.exists(gdp_growth_path)) {
  stop(gdp_growth_path, " not found - run Rscript R/01_fetch_data.R first.",
       call. = FALSE)
}
gdp_growth <- readRDS(gdp_growth_path)

dir.create("data/raw",       showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# --- Hours worked: ABS 5206.0 Table 1 (Key Aggregates) -----------------------
raw_t1 <- with_retry(function() read_abs("5206.0", tables = 1),
                     "read_abs 5206.0 table 1")

hours <- raw_t1 |>
  filter(
    series_type == "Seasonally Adjusted",
    series == "Hours worked: Index ;"
  ) |>
  select(date, hours_idx = value) |>
  arrange(date)

stopifnot(nrow(hours) > 0)

hours_growth <- hours |>
  mutate(g_U = annualised_qoq(hours_idx)) |>
  select(date, g_U) |>
  filter(!is.na(g_U))

# Merge with gdp_growth on date; m8 requires all four series at every t.
gdp_growth_4 <- gdp_growth |>
  inner_join(hours_growth, by = "date")

saveRDS(gdp_growth_4, file = "data/processed/gdp_growth_with_hours.rds")

# --- Monthly employment: ABS 6202.0 Table 1 (Labour Force Survey) ------------
# Used by m11 (mixed-frequency Banbura-Modugno) as a monthly activity indicator.
lfs_files <- with_retry(function() show_available_files("labour-force-australia"),
                        "show_available_files labour-force-australia")
lfs_url   <- lfs_files$url[lfs_files$file == "62020001.xlsx"]
stopifnot(length(lfs_url) == 1)
lfs_raw <- with_retry(function() read_abs_url(lfs_url), "read_abs_url 62020001.xlsx")

employment <- lfs_raw |>
  filter(
    series_type == "Seasonally Adjusted",
    series == "Employed total ;  Persons ;"
  ) |>
  select(date, employed = value) |>
  arrange(date)
stopifnot(nrow(employment) > 100)

# Month-on-month annualised growth: 1200 x dlog
emp_growth <- employment |>
  mutate(g_emp = 1200 * (log(employed) - log(lag(employed)))) |>
  filter(!is.na(g_emp)) |>
  select(date, g_emp)

saveRDS(employment, file = "data/processed/employment_monthly.rds")
saveRDS(emp_growth, file = "data/processed/employment_growth_monthly.rds")

cat(sprintf(
  "Hours-worked series available for T = %d quarters (m8 uses this subset).\n",
  nrow(gdp_growth_4)
))
cat(sprintf(
  "Monthly employment growth: %d months from %s to %s (m11 uses this).\n",
  nrow(emp_growth), min(emp_growth$date), max(emp_growth$date)
))
