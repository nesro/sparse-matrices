#!/bin/bash
#
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2013-2014
# https://github.com/nesro/sparse-matrices
#
# This is a helper file for running jobs on our school server STAR.
# Please do NOT to try hack this in any way.

set -x

./tests/graphs/mul_speed/run.sh
exit 0

make
time ./main -f csr -a <( gzip -cd ../test_matrices/generated/test1.mtx.gz) -v
time ./main -f coo -a <( gzip -cd ../test_matrices/generated/test1.mtx.gz) -v

exit 0



matrix_a=../test_matrices/test8/test8a.mtx.gz
matrix_b=../test_matrices/test8/test8b.mtx.gz

echo "--csr--"
time ./main -f csr -a <( gzip -cd $matrix_a ) -b <( gzip -cd $matrix_b ) -v

for i in 2; do
for j in 256; do

echo "i=$i j=$j"
make KAT_N=$i
time ./main -f kat -s $j -a <( gzip -cd $matrix_a ) -b <( gzip -cd $matrix_b )  -v

echo " ... "

done
done
