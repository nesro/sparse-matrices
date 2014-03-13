#!/bin/bash

i=$(cat ./tests/eia/cachegrind/last_O0NO.txt)
echo $((i + 512)) > ./tests/eia/cachegrind/last_O0NO.txt

make clean 2>&1 >/dev/null
DEBUG=1 make 2>&1 >/dev/null
valgrind --tool=cachegrind ./main cn $i a 2>&1 >/dev/null
#Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
totals=$(callgrind_annotate --auto=yes cachegrind* | grep "PROGRAM TOTALS" | tr -d ',')
echo "$i $totals" > ./tests/eia/cachegrind/O0NO_i_${i}_`date +%s`.txt
