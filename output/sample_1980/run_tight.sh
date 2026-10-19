#!/bin/bash
# Driver for the tight-priors variant of m13 (stan/m13_tight_priors.stan; NOTES.md,
# "Tighter priors for m13", 2026-10-19). Sequential, resumable, one Stan job at a time:
#   nohup bash output/sample_1980/run_tight.sh >> output/sample_1980/logs/driver_tight.log 2>&1 &
# Each step is skipped if its output already exists. Production files are not touched.
cd /Users/davidstephan/Documents/latentgdp/LatentGDP || exit 1
L=output/sample_1980/logs
stamp() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*"; }
stamp "DRIVER START"
for s in 1 2; do
  if [ "$s" = 1 ]; then f=output/sample_1980/fit_m13_tight_gate.rds; else f=output/sample_1980/fit_m13_tight_seed2_gate.rds; fi
  if [ -e "$f" ]; then stamp "skip seed $s ($f exists)"; continue; fi
  stamp "start m13_tight seed $s"
  LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980 \
    Rscript R/fit_m13.R --variant=tight --seed=$s > $L/m13_tight_seed$s.log 2>&1
  stamp "end m13_tight seed $s (exit $?)"
done
if [ -e output/sample_rdp/fit_m13_tight.rds ]; then stamp "skip rdp (exists)"; else
  stamp "start m13_tight rdp sample"
  Rscript paper/fit_rdp_sample.R m13_tight > output/sample_rdp/logs/fit_rdp_m13_tight.log 2>&1
  stamp "end m13_tight rdp sample (exit $?)"
fi
stamp "DRIVER DONE"
