#!/bin/bash

server=$1
action=$2
format=$3
matrix_size=$4
density=$5
random_seed=$6

gp_mul_time
gp_load_a_time_fi

action_list = compute compute_all generate_graph

function usage {
	echo "Usage $0 [action] [format] [matrix_size] [density] [random_seed]
}

for (( i = 128; i <= 8196; i*=2)); do
	
done
