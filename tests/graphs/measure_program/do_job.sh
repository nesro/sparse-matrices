#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

#-------------------------------------------------------------------------------
# read input
#-------------------------------------------------------------------------------

matrix_dir=${1:-0}
matrix_name=${2:-0}
mul_vector=${3:-0}
format=${4:-0}
block_size=${5:-0}
kat_n=${6:-0}

echo "matrix_dir $matrix_dir"
echo "matrix_name $matrix_name"
echo "mul_vector $mul_vector"
echo "format $format"
echo "block_size $block_size"
echo "kat_n $kat_n"

#-------------------------------------------------------------------------------
# check input
#-------------------------------------------------------------------------------

if [[ "$matrix_dir" == "0" ]]; then
	echo "matrix a must be set"
	exit 1
fi

if [[ "$matrix_dir" == "0" ]]; then
	echo "matrix a must be set"
	exit 1
fi

if [[ "$format" == "0" ]]; then
	echo "format must be set"
	exit 1
fi

#-------------------------------------------------------------------------------
# find matices
#-------------------------------------------------------------------------------

matrix_a=${matrix_dir}/${matrix_name}.mtx.gz

if [[ ! -f $matrix_a ]];
	echo "matrix_a ( $matrix_a ) not found!"
	exit 1
fi

if [[ "$mul_vector" == "1" ]]; then
	matrix_b=${matrix_dir}/vector_${matrix_name}.mtx.gz
	if [[ ! -f $matrix_b ]];
		echo "vector matrix_b ( $matrix_b ) not found!"
		exit 1
	fi
	outfile=out_mvm_${matrix_name}_${format}_${block_size}_${kat_n}_.txt
else
	outfile=out_mmm_${matrix_name}_${format}_${block_size}_${kat_n}_.txt
fi

#-------------------------------------------------------------------------------
# run the computation
#-------------------------------------------------------------------------------

set -x

make PRECISION=2 KAT_N=$kat_n

if (( $matrix_b == 0 )); then
	# matrix * matrix
	time ./main -f $format $block_size_param \
		-a <( gzip -cd $matrix_a ) -v >$outfile
else
	# matrix * vector
	time ./main -f $format $block_size_param -V \
		-a <( gzip -cd $matrix_a ) -b <( gzip -cd $matrix_b ) -v >$outfile
fi

#-------------------------------------------------------------------------------
