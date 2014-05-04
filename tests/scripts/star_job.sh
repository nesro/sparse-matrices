#!/bin/bash
#
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
# https://github.com/nesro/sparse-matrices
#
# This is a helper file for running jobs on our school server STAR.
# Please do NOT to try hack this in any way.

set -x

matrix=../test_matrices/test6.mtx.gz

echo "--csr--"
time ./main -f csr -a <( gzip -cd $matrix ) -v

for i in 2 4; do
for j in 256 512 1024; do

echo "i=$i j=$j"
make KAT_N=$i
time ./main -f kat -s $j -a <( gzip -cd $matrix ) -v

echo " ... "

done
done
