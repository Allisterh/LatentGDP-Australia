# Production data fetch: ABS 5206.0 Table 24 only (GDP-E/I/P chain volume
# measures). This is all the production model (m9) and the weekly CI job need.
# The research-only series (5206.0 Table 1 hours worked for m8/m10, and the
# 6202.0 monthly employment for m11) are fetched by R/01b_fetch_research_data.R.
#
# Writes:
#   data/processed/gdp_eip_levels.rds  - long levels: date, approach, value
#   data/processed/gdp_growth.rds      - wide growth: date, g_E, g_I, g_P

suppressPackageStartupMessages({
  library(readabs)
  library(dplyr)
  library(tidyr)
})

source("R/utils_abs.R")
neutralise_abs_connection_check()

dir.create("data/raw",       showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# ABS Cat. 5206.0 Table 24 = Selected Analytical Series (GDP-E/I/P CVMs).
raw <- with_retry(function() read_abs("5206.0", tables = 24),
                  "read_abs 5206.0 table 24")

eip_pattern <- paste0(
  "^Gross domestic product - (Expenditure|Income|Production)",
  " based: Chain volume measures ;$"
)

gdp_eip <- raw |>
  filter(
    series_type == "Seasonally Adjusted",
    grepl(eip_pattern, series)
  ) |>
  mutate(
    approach = case_when(
      grepl("Expenditure", series) ~ "E",
      grepl("Income",      series) ~ "I",
      grepl("Production",  series) ~ "P"
    )
  ) |>
  select(date, approach, value) |>
  arrange(date, approach)

# Quarterly log-difference, annualised (x400) - matches paper units.
gdp_growth <- gdp_eip |>
  group_by(approach) |>
  arrange(date, .by_group = TRUE) |>
  mutate(g = annualised_qoq(value)) |>
  ungroup() |>
  select(date, approach, g) |>
  pivot_wider(names_from = approach, values_from = g, names_prefix = "g_") |>
  drop_na()

stopifnot(all(c("g_E", "g_I", "g_P") %in% names(gdp_growth)))

saveRDS(gdp_eip,    file = "data/processed/gdp_eip_levels.rds")
saveRDS(gdp_growth, file = "data/processed/gdp_growth.rds")

cat(sprintf(
  "Saved %d quarters of levels from %s to %s; growth T = %d quarters from %s to %s.\n",
  nrow(gdp_eip) / 3, min(gdp_eip$date), max(gdp_eip$date),
  nrow(gdp_growth), min(gdp_growth$date), max(gdp_growth$date)
))
