#!/bin/bash
# Driver for the paper's new analyses (sequential, resumable). Run from the repo root:
#   nohup bash output/sample_rdp/run_new_analyses.sh >> output/sample_rdp/logs/driver.log 2>&1 &
cd "$(dirname "$0")/../.."
echo "[$(date '+%F %T')] DRIVER START"
Rscript paper/fit_rdp_sample.R m2 m13 m5 m4 m9 m12 > output/sample_rdp/logs/fit_rdp.log 2>&1
echo "[$(date '+%F %T')] rdp fits exit $?"
if [ ! -f output/sample_1980/fit_m13_omegag1.rds ]; then
  LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980 Rscript R/fit_m13.R --omega_g_sd=1 \
    --out=output/sample_1980/fit_m13_omegag1.rds > output/sample_1980/logs/m13_omegag1.log 2>&1
  echo "[$(date '+%F %T')] omega_g sensitivity fit exit $?"
fi
echo "[$(date '+%F %T')] DRIVER DONE"
