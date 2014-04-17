#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

source ./tests/scripts/utils.sh

prog=./tests/bin/mulres_distance
mtx_dir=./tests/matrices/generated
mm_dir=./tests/matrices/matrix_market
gp=./tests/graphs/strassen_stability/gp

#echo "orsirr_1">$gp.txt
#for ((i=1;i<=512;i*=2)); do
#	echo $i
#	echo "$($prog $mm_dir/orsirr_1.mtx $mm_dir/orsirr_1.mtx $i) $i">>$gp.txt
#done

gnuplot_wrapper $gp.txt $gp.pdf "práh velikosti podmatice" "suma rozdílu součtů" 2
