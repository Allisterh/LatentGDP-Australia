#!/bin/bash
# Added after the main driver: COVID particle checks for m9 and m4, whose 1980-sample
# filter passes disagree by 1.5-4 nats (the 1980 fits treat 2020 Q2 as common measurement
# error, as m12's do). Resumable like run_lfo_1980.sh.
cd /Users/davidstephan/Documents/latentgdp/LatentGDP || exit 1
export LATENTGDP_START=1980-03-01 LATENTGDP_OUT=output/sample_1980
LOG=output/sample_1980/lfo
step() { local name=$1; shift; echo "$(date '+%F %T') START $name: Rscript R/03_compare.R $*"
  Rscript R/03_compare.R "$@" --cores=4 >> "$LOG/run_$name.log" 2>&1; echo "$(date '+%F %T') END $name exit $?"; }
step pfcheck_covid_m9m4 pfcheck m9 m4
step pfref_m9m4 pfref m9 m4
step summary summary
echo "$(date '+%F %T') EXTRA DONE"
