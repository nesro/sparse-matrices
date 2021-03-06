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
0)
	matrix_dir=${2:-0}
	matrices=${3:-0}
	;;
1)
	matrix_dir="../test_matrices/generated"
	matrices="test1 test2 test3 test4"
	;;
2)
	matrix_dir="../test_matrices/final"
	matrices="fp EX6 gupta3 appu human_gene2 c8_mat11_I heart1 exdata_1"
	;;
3)
	matrix_dir="../test_matrices/dense"
	matrices="d256 d512 d1024 d2048"
	;;
*)
	echo "a job type must be set"
	exit 1
	;;
esac

#------------------------------------------------------------------------------

# mtx_dir mtx_name mul_vector format block_size kat_n
script=./tests/scripts/star.sh

for mul_vector in 0 1; do
	for matrix in $matrices; do
		for f in coo csr; do
		
			# Let's skip SpMMM in the COO format
			if [[ "$f" == "coo" && "$mul_vector" == "0" ]]; then
				continue;
			fi
		
			$script $matrix_dir $matrix $mul_vector $f 0 2
			wait_for_slot
		done
		for bs in 2 4 8 16 32 64; do
			for f in bsr; do
				$script $matrix_dir $matrix $mul_vector $f $bs 2
				wait_for_slot
			done
			for f in kat; do
				for k in 2 4 8; do
					$script $matrix_dir $matrix $mul_vector $f $bs $k
					wait_for_slot
				done
			done
		done
	done
done

echo "create jobs completed"
