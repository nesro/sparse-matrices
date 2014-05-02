#!/bin/bash
#
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
# https://github.com/nesro/sparse-matrices
#
# This is a helper file for running jobs on our school server STAR.
# Please do NOT to try hack this in any way.

make KAT_N=4 tests
#time ./main -f csr -a <(gzip -cd ./matrices/cage12.mtx.gz) -b <(gzip -cd ./tests/matrices/vector_cage12_130228.mtx.gz) -V -v -o ./out.mtx
time ./tests/graphs/format_size/run.sh