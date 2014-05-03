#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

set -x

source ./tests/scripts/utils.sh

make DOUBLE_PRECISION=1

# time ./main -f kat -s 2 \
# -a <(gzip -cd ../big_matrices/ldoor.mtx.gz) \
# -b <(gzip -cd ../big_vectors/vector_ldoor_952203.mtx.gz) \
# -V -v -o ./resvec_ldoor.mtx

if false; then
	./main -f bsr -s 64 \
	-a <(gzip -cd ../big_matrices/cage12.mtx.gz) \
	-b <(gzip -cd ../big_vectors/vector_cage12_130228.mtx.gz) \
	-V -v -o ./resvec_ldoor.mtx
fi

time for matrix in $big_list; do

	for format in \
		"coo" \
		"csr" \
		"bsr -s 64" \
		"bsr -s 128" \
		"bsr -s 256" \
		"kat -s 64" \
		"kat -s 128" \
		"kat -s 256"; do
		
		echo -n "\"$matrix\" \"$format\"" >> gp_$$.txt
		
		time ./main -f ${format} \
			-a <(gzip -cd ${big_dir_mat}/${matrix}.mtx.gz) \
			-b <(gzip -cd ${big_dir_vec}/vector_${matrix}_*.mtx.gz) \
			-V -v -o ./resvec_${matrix}_$(echo ${format} | tr ' ' '_').mtx | grep "a_size" | cut -d' ' -f2 >> gp_$$.txt
		
		echo "CHECK"
		for check in $(ls ./resvec_${matrix}_*); do
			diff -q $check ./resvec_${matrix}_$(echo ${format} | tr ' ' '_').mtx
		done
		
	done
done
