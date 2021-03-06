#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

set -x

source ./tests/scripts/utils.sh


make PRECISION=2

echo -n "format "> gp_$$_res.txt
for matrix in $generated_list; do echo "$matrix " | tr -d '\n' >>gp_$$_res.txt; done
echo " ">>gp_$$_res.txt

katon=0
bsron=1


# begin formats
if (( $bsron == 1 )); then
for format in \
		"csr" \
		"bsr -s 16" \
		"bsr -s 32" \
		"bsr -s 64" \
		"bsr -s 128" \
		"bsr -s 256"; do
		
		echo -n "\"${format}\" " >gp_$$_$(echo $format | tr ' ' '_').txt
done
fi

if (( $katon == 1 )); then
for i in 2 4; do
	for format in \
				"kat -s 16" \
				"kat -s 32" \
				"kat -s 64" \
				"kat -s 128" \
				"kat -s 256"; do
		ktmp="k=$i ${format}"
		echo -n "\"$ktmp\" " >>gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt
	done
done 
fi

if (( $bsron == 1 )); then
time for matrix in $generated_list; do

	for format in \
		"csr" \
		"bsr -s 16" \
		"bsr -s 32" \
		"bsr -s 64" \
		"bsr -s 128" \
		"bsr -s 256"; do

		time ./main -f ${format} \
			-a <(gzip -cd ${generated_dir}/${matrix}.mtx.gz) \
			-v | grep "a_size" | cut -d' ' -f2 | tr -d '\n' >>gp_$$_$(echo $format | tr ' ' '_').txt
		
		echo -n " ">>gp_$$_$(echo $format | tr ' ' '_').txt
	done
done
fi

if (( $katon == 1 )); then
for i in 2 4; do
	make PRECISION=2 KAT_N=$i
	time for matrix in $generated_list; do
		for format in \
				"kat -s 16" \
				"kat -s 32" \
				"kat -s 64" \
				"kat -s 128" \
				"kat -s 256"; do
			
			ktmp="k=$i ${format}"
					
			time ./main -f ${format} \
			-a <(gzip -cd ${generated_dir}/${matrix}.mtx.gz) \
			-v | grep "a_size" | cut -d' ' -f2 | tr -d '\n' >>gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt
			echo -n " " >>gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt
		done
	done
done
fi

# end formats
if (( $bsron == 1 )); then
for format in \
		"csr" \
		"bsr -s 16" \
		"bsr -s 32" \
		"bsr -s 64" \
		"bsr -s 128" \
		"bsr -s 256"; do
		
		echo " " >>gp_$$_$(echo $format | tr ' ' '_').txt
		cat gp_$$_$(echo $format | tr ' ' '_').txt	>>gp_$$_res.txt
done
fi

if (( $katon == 1 )); then
for i in 2 4; do
	for format in \
				"kat -s 16" \
				"kat -s 32" \
				"kat -s 64" \
				"kat -s 128" \
				"kat -s 256"; do
		ktmp="k=$i ${format}"
		echo " " >>gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt
		cat gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt >>gp_$$_res.txt
	done
done 
fi

