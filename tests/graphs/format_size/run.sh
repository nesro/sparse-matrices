#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

set -x

source ./tests/scripts/utils.sh

for matrix in $big_list; do

	for format in \
		"coo" \
		"csr" \
		"bsr -s 64" \
		"bsr -s 128" \
		"bsr -s 256" \
		"kat -s 64" \
		"kat -s 128" \
		"kat -s 256"; do
		
		time ./main -f ${format} \
			-a <(gzip -cd ${big_dir_mat}/${matrix}.mtx.gz) \
			-b <(gzip -cd ${big_dir_vec}/vector_${matrix}_*.mtx.gz) \
			-V -v -o ./resvec_${matrix}.mtx >>log_$$.txt
		
	done
done
