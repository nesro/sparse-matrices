#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

if (( $# != 1 )); then
	echo "Usage: $0 n"
	exit 1
fi

for ((i=0;i<$1;i++)); do
	for ((j=0;j<$1;j++)); do
		for ((k=0;k<$1;k++)); do
			echo "IN_KAT_UNROLL_MUL(A,B,C,$i,$k,$k,$j); \\"
		done
	done
done