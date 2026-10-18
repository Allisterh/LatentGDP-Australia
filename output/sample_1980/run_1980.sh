#!/bin/bash
# Driver for the 1980 Q1 sample (RBA RDP 2014-12 start): sequential, resumable.
#   nohup bash output/sample_1980/run_1980.sh > output/sample_1980/logs/driver.log 2>&1 &
# Each step is skipped if its output already exists. One heavy job at a time.
cd /Users/davidstephan/Documents/latentgdp/LatentGDP || exit 1
export LATENTGDP_START=1980-03-01
export LATENTGDP_OUT=output/sample_1980
O=$LATENTGDP_OUT; L=$O/logs
stamp() { echo "$(date '+%Y-%m-%d %H:%M:%S') $*"; }
step() {  # step NAME DONE_FILE CMD...
  local name=$1 done=$2; shift 2
  if [ -e "$done" ]; then stamp "skip $name ($done exists)"; return 0; fi
  stamp "start $name"
  "$@" > "$L/$name.log" 2>&1; local rc=$?
  stamp "end $name (exit $rc)"
  return 0
}
has_row() { [ -f $O/tables/fit_diagnostics.csv ] && grep -q "^\"$1\"," $O/tables/fit_diagnostics.csv; }
fit02() {  # fit one model through R/02_fit_models.R unless fit + row exist
  local m=$1
  if [ -e $O/fit_$m.rds ] && has_row $m; then stamp "skip fit_$m (fit + row exist)"; return 0; fi
  stamp "start fit_$m"; Rscript R/02_fit_models.R $m > $L/fit_$m.log 2>&1; stamp "end fit_$m (exit $?)"
}
for m in m1 m2 m5 m6 m7 m8 m3; do fit02 $m; done
step m13_seed1 $O/fit_m13_gate.rds Rscript R/fit_m13.R
step m13_seed2 $O/fit_m13_seed2_gate.rds Rscript R/fit_m13.R --seed=2
if [ -e $O/fit_m13.rds ] && ! has_row m13; then
  stamp "diag m13"; Rscript R/02_fit_models.R --diag-only m13 > $L/diag_m13.log 2>&1; stamp "end diag m13 (exit $?)"; fi
for m in m4 m9 m10 m12; do fit02 $m; done
step fit_m11 $O/fit_m11.rds Rscript R/fit_m11.R
if [ -e $O/fit_m11.rds ] && ! has_row m11; then
  stamp "diag m11"; Rscript R/02_fit_models.R --diag-only m11 > $L/diag_m11.log 2>&1; stamp "end diag m11 (exit $?)"; fi
stamp "vardecomp"; Rscript R/02_fit_models.R --vardecomp > $L/vardecomp.log 2>&1; stamp "end vardecomp (exit $?)"
stamp "figures"; Rscript R/04_figures.R > $L/figures.log 2>&1; stamp "end figures (exit $?)"
stamp "DRIVER DONE"
