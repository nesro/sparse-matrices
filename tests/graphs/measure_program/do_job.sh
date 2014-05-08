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
kat_n=${6:-2}

use_gzipped="1"

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

matrix_a=${matrix_dir}/${matrix_name}.mtx
matrix_b=0

#if [[ ! -f $matrix_a ]]; then
#	echo "matrix_a ( $matrix_a ) not found!"
#	exit 1
#fi

safe_matrix_name=$(echo $matrix_name | tr '_' '-')

if [[ "$mul_vector" == "1" ]]; then
	matrix_b=${matrix_dir}/vector_${matrix_name}.mtx
#	if [[ ! -f $matrix_b ]]; then
#		echo "vector matrix_b ( $matrix_b ) not found!"
#		exit 1
#	fi
	outfile=out_mvm_${safe_matrix_name}_${format}_${block_size}_${kat_n}_.txt
else
	outfile=out_mmm_${safe_matrix_name}_${format}_${block_size}_${kat_n}_.txt
fi

if [[ $block_size == "0" ]]; then
	block_size_param=""
else
	block_size_param="-s $block_size"
fi

#-------------------------------------------------------------------------------
# run the computation
#-------------------------------------------------------------------------------

set -x

make PRECISION=2 KAT_N=$kat_n
if  [[ "$use_gzipped" == "1" ]]; then

	echo "using gzipped matrix_b = ${matrix_b}"

	if [[ "$matrix_b" == "0" ]]; then
		# matrix * matrix
		time ./main -f $format $block_size_param \
			-a <( gzip -cd ${matrix_a}.gz ) -v 2>&1 >$outfile
	else
		# matrix * vector
		time ./main -f $format $block_size_param -r 100 -V \
			-a <( gzip -cd ${matrix_a}.gz ) -b <( gzip -cd ${matrix_b}.gz ) -v 2>&1 >$outfile
	fi
else # use_gzipped

	echo "using unzipped matrix_b = ${matrix_b}"

	if [[ "$matrix_b" == "0" ]]; then
		# matrix * matrix
		time ./main -f $format $block_size_param \
			-a $matrix_a -v 2>&1 >$outfile
	else
		# matrix * vector
		time ./main -f $format $block_size_param -r 100 -V \
			-a $matrix_a -b $matrix_b -v 2>&1 >$outfile
	fi
fi # use_gzipped

#-------------------------------------------------------------------------------
