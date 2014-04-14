#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

# ./tests/matrix_generator.c
generator=./tests/matrix_generator

if [[ ! -x $generator ]]; then
	echo "Matrix generator was not found at \"$generator\"." \
		"Maybe you want to run make?"
	exit 1
fi

function generate {
	file="./tests/matrices/generated/$1.mtx"
	shift
	
	#if [[ ! -e $file ]]; then	
		$generator -o $file $@	
	#fi
}

mkdir -p ./tests/matrices/generated

generate "f64x64_01" -n 64 -S 100000.009
generate "f64x64_02" -n 64 -S 1

generate "16x16_01" -n 16 -S -200
generate "16x16_02" -n 16 -S -201

generate "32x32_01" -n 32 -S -200
generate "32x32_02" -n 32 -S -201

generate "64x64_01" -n 64 -S -200
generate "64x64_02" -n 64 -S -201

generate "128x128_01" -n 128 -S -200
generate "128x128_02" -n 128 -S -201

generate "256x256_01" -n 256 -S -200
generate "256x256_02" -n 256 -S -201

generate "512x512_01" -n 512 -S -200
generate "512x512_02" -n 512 -S -201
