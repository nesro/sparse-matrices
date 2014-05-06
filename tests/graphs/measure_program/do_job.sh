#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

matrix_a=${1:-0}
matrix_b=${2:-0} # 0 if not set
format=${3:-0}
block_size=${4:-0}
kat_n=${5:-2}

echo "matrix_a $matrix_a"
echo "matrix_b $matrix_b"
echo "format $format"
echo "block_size $block_size"
echo "kat_n $kat_n"

if [[ "$matrix_a" == "0" ]]; then
	echo "matrix a must be set"
	exit 1
fi

if [[ "$format" == "0" ]]; then
	echo "format must be set"
	exit 1
fi

if [[ "$block_size" == "0" ]]; then
	block_size_param=""
else
	block_size_param="-s $block_size"
fi

matrix_a_file=$(echo $matrix_a | rev |  cut -d'/' -f1 | rev | cut -d'.' -f1)

if [[ "$matrix_b" == "0" ]]; then
	matrix_b_file="-"
else
	matrix_b_file=$(echo $matrix_b | rev |  cut -d'/' -f1 | rev | cut -d'.' -f1)
fi

outfile=out_${matrix_a_file}_${matrix_b_file}_${format}_${block_size}_${kat_n}.txt

make PRECISION=2 KAT_N=$kat_n

echo $outfile

set -x

if (( $matrix_b == 0 )); then
	time ./main -f $format $block_size_param \
		-a <( gzip -cd $matrix_a ) -v >$outfile
else
	time ./main -f $format $block_size_param -V \
		-a <( gzip -cd $matrix_a ) -b <( gzip -cd $matrix_b ) -v >$outfile
fi

