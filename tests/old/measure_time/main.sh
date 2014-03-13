#!/bin/bash

server=$1
lo=$2
hi=$3

formats=${4:"csr quadtree"}

binary=./main

gp_mul_time
gp_load_a_time_fi

action_list = compute compute_all generate_graph

function usage {
	echo "Usage $0 [action] [format] [matrix_size] [density] [random_seed]
}

for (( i = lo; i <= hi; i*=2 )); do

done

for (( i = lo; i <= hi; i*=2 )); do

	# do matrix a and b exists?

	for format in $formats; do
		res=$($binary -f $format -a $matrix_a -b $matrix_b -o none)
		mul=$(grep multiplication <($res))
		load=$(grep load_a <($res))
	done
done
