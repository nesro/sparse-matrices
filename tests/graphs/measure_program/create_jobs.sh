#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

source ./tests/scripts/utils.sh

set -x

jobs_type=${1:-0}

matrix_dir=0
matrices=0
mul_vector=0

case $jobs_type in
1)
	matrix_dir="../test_matrices/generated"
	matrices="test1 test2 test3 test4"
	;;
2)
	matrix_dir="../test_matrices/generated"
	matrices="test1 test2 test3 test4"
	mul_vector=1	
	;;
*)
	echo "a job type must be set"
	exit 1
	;;
esac

#------------------------------------------------------------------------------

# mtx_dir mtx_name mul_vector format block_size kat_n
script=./tests/scripts/star.sh

for matrix in $matrices; do
	for f in csr; do
		$script $matrix_dir $matrix $mul_vector $f 0 2
		wait_for_slot
	done
	for bs in 16 32 64 128 256 512; do
		for f in bsr; do
			$script $matrix_dir $matrix $mul_vector $f $bs 2
			wait_for_slot
		done
		for f in kat; do
			for k in 2 4; do
				$script $matrix_dir $matrix $mul_vector $f $bs $k
				wait_for_slot
			done
		done
	done
done

echo "create jobs completed"
