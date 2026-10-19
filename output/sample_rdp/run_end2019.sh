#!/bin/bash
# Fits ending in 2019 Q4 (after run_new_analyses.sh has finished). Run from the repo root:
#   nohup bash output/sample_rdp/run_end2019.sh >> output/sample_rdp/logs/driver.log 2>&1 &
cd "$(dirname "$0")/../.."
until grep -q "DRIVER DONE" output/sample_rdp/logs/driver.log; do sleep 30; done
echo "[$(date '+%F %T')] END2019 START"
Rscript paper/fit_rdp_sample.R --end=2019-12-01 --out=output/sample_rdp/end2019 m2 m13 > output/sample_rdp/logs/fit_end2019.log 2>&1
echo "[$(date '+%F %T')] END2019 exit $?"
