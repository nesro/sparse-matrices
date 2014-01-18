#!/bin/bash

LIMIT=8192

if [[ $# -eq 1 ]]; then
  LIMIT=$1
fi

echo "- loop_unrolling gcc_opt loop_unrolling+gcc_opt" > results

for i in $(seq 512 512 $LIMIT); do
  echo "(./tests/csr.sh MVM TEST $i | grep time | cut -d= -f2 | cut -d' ' -f2 | tr '\n' ' ' && echo $i) >> results"
  (./tests/csr.sh MVM TEST $i | grep time | cut -d= -f2 | cut -d' ' -f2 | tr '\n' ' ' && echo $i) >> results
done

cat gnuplot.conf | gnuplot
