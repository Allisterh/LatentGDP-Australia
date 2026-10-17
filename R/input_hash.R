# Input fingerprint for the production pipeline.
#
# SHA-256 over (a) a canonical text serialisation of the E/I/P growth data the
# production model is fitted to and (b) the bytes of every file that defines
# the production pipeline (Stan model, R scripts, CI workflow). If none of
# these change, a re-run would reproduce the published estimate, so the weekly
# CI job compares this hash with the "input_hash" field of the committed
# docs/diagnostics.json and skips fit / build / commit / deploy when equal.
#
# Usage:
#   Rscript R/input_hash.R            # prints the hash (one line) to stdout
#   source("R/input_hash.R"); compute_input_hash()
#
# Requires the 'openssl' package (a dependency of readabs via httr).

INPUT_HASH_FILES <- c(
  "stan/m12_struct_sv_serial.stan",
  "R/utils_abs.R",
  "R/01_fetch_data.R",
  "R/input_hash.R",
  "R/fit_production.R",
  "R/output_estimate.R",
  "R/build_site.R",
  ".github/workflows/recompute.yml"
)

# Canonical CSV text of the model inputs: fixed column order, ISO dates,
# fixed 6-decimal formatting (annualised % units), LF line endings. Six
# decimals is far below any ABS revision yet robust to last-bit libm
# differences between platforms.
canonical_gdp_text <- function(gdp_growth) {
  stopifnot(all(c("date", "g_E", "g_I", "g_P") %in% names(gdp_growth)))
  rows <- sprintf("%s,%.6f,%.6f,%.6f",
                  format(as.Date(gdp_growth$date), "%Y-%m-%d"),
                  gdp_growth$g_E, gdp_growth$g_I, gdp_growth$g_P)
  paste(c("date,g_E,g_I,g_P", rows), collapse = "\n")
}

sha256_text <- function(x) {
  as.vector(unclass(as.character(openssl::sha256(charToRaw(enc2utf8(x))))))
}

sha256_file_lf <- function(path) {
  if (!file.exists(path)) stop("input_hash: missing file ", path, call. = FALSE)
  raw <- readBin(path, "raw", n = file.info(path)$size)
  txt <- rawToChar(raw)
  txt <- gsub("\r\n", "\n", txt, fixed = TRUE)   # normalise CRLF checkouts
  sha256_text(txt)
}

input_hash_components <- function(data_path = "data/processed/gdp_growth.rds",
                                  files = INPUT_HASH_FILES) {
  if (!file.exists(data_path)) {
    stop("input_hash: missing ", data_path, " - run R/01_fetch_data.R first.",
         call. = FALSE)
  }
  gdp_growth <- readRDS(data_path)
  c(
    stats::setNames(sha256_text(canonical_gdp_text(gdp_growth)), "data:gdp_growth"),
    stats::setNames(vapply(files, sha256_file_lf, character(1)),
                    paste0("file:", files))
  )
}

compute_input_hash <- function(...) {
  comp <- input_hash_components(...)
  sha256_text(paste(names(comp), comp, sep = " ", collapse = "\n"))
}

if (sys.nframe() == 0L) {
  cat(compute_input_hash(), "\n", sep = "")
}
