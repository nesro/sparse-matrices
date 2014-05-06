#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

matrix_set=${1:-0}
vector_set=${2:-0}

set -x

sleep_time=45

if (( $matrix_set == 0 )); then
	matrix_set="../test_matrices/generated/test1.mtx.gz \
../test_matrices/generated/test2.mtx.gz \
../test_matrices/generated/test3.mtx.gz \
../test_matrices/generated/test4.mtx.gz"
fi

#if (( $vector_set != 0 )); then
	#todo
#fi 

# matA matB format block_size kat_n
# script=./tests/graphs/measure_program/do_job.sh
script=./tests/scripts/star.sh

for matrix in $matrix_set; do
	for vector in $vector_set; do
		for f in csr; do
			$script $matrix $vector $f 0 2
		done
		for bs in 16 32 64 128 256 512; do
			for f in bsr; do
				$script $matrix $vector $f $bs 2
				sleep $sleep_time
			done
			for f in kat; do
				for k in 2 4; do
					$script $matrix $vector $f $bs $k
					sleep $sleep_time
				done
			done
		done
	done
done

echo "create jobs completed"
