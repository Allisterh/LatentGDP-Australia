#!/bin/bash
# PSIS-LFO-CV on the 1980 Q1 sample (model-comparison role, 2026-10-18). Sequential and
# resumable: every R/03_compare.R step resumes from its checkpoints in output/sample_1980/lfo,
# so the whole script can simply be re-run after an interruption.
#   nohup bash output/sample_1980/lfo/run_lfo_1980.sh >> output/sample_1980/lfo/driver.log 2>&1 &
cd /Users/davidstephan/Documents/latentgdp/LatentGDP || exit 1
export LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980
LOG=output/sample_1980/lfo
step() {   # step <name> <args...>
  local name=$1; shift
  echo "$(date '+%F %T') START $name: Rscript R/03_compare.R $*"
  Rscript R/03_compare.R "$@" --cores=4 >> "$LOG/run_$name.log" 2>&1
  local rc=$?
  echo "$(date '+%F %T') END $name exit $rc"
}
# linear-Gaussian models (exact Kalman filter), then the SV models (main pass + 3 re-runs each)
for m in m1 m2 m3 m5 m6 m7 m8 m11; do step $m $m; done
for m in m13 m12 m9; do step $m $m; done
# particle-filter bias checks: COVID targets, extra 2020Q2 references, 2015-19 targets
step pfcheck_covid pfcheck m13 m12
step pfref pfref m13 m12
step pfcheck_pre pfcheck m13 m12 --precovid
for m in m4 m10; do step $m $m; done
# seed-selection (retry) sensitivity: every model; models without retried fits return at once
step retrycheck retrycheck m1 m2 m3 m4 m5 m6 m7 m8 m9 m10 m11 m12 m13
step summary summary
echo "$(date '+%F %T') DRIVER DONE"
