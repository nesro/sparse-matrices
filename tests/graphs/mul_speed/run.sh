#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

set -x

source ./tests/scripts/utils.sh

make DOUBLE_PRECISION=1



echo -n "format "> gp_$$_res.txt
for matrix in $big_list; do echo "$matrix " | tr -d '\n' >>gp_$$_res.txt; done
echo " ">>gp_$$_res.txt

# begin formats
for format in \
		"coo" \
		"csr" \
		"bsr -s 4" \
		"bsr -s 8" \
		"bsr -s 16"; do
		
		echo -n "\"${format}\" " >gp_$$_$(echo $format | tr ' ' '_').txt
done
for i in 2 4; do
	for format in \
				"kat -s 4" \
				"kat -s 8" \
				"kat -s 16"; do
		ktmp="k=$i ${format}"
		echo -n "\"$ktmp\" " >>gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt
	done
done 



time for matrix in $big_list; do

	for format in \
		"coo" \
		"csr" \
		"bsr -s 4" \
		"bsr -s 8" \
		"bsr -s 16"; do

		time ./main -f ${format} \
			-a <(gzip -cd ${big_dir_mat}/${matrix}.mtx.gz) \
			-v | grep "time_mul" | cut -d' ' -f2 | tr -d '\n' >>gp_$$_$(echo $format | tr ' ' '_').txt
		
		echo -n " ">>gp_$$_$(echo $format | tr ' ' '_').txt
	done
done

for i in 2 4; do
	make DOUBLE_PRECISION=1 KAT_N=$i
	time for matrix in $big_list; do
		for format in \
				"kat -s 4" \
				"kat -s 8" \
				"kat -s 16"; do
			
			ktmp="k=$i ${format}"
					
			time ./main -f ${format} \
			-a <(gzip -cd ${big_dir_mat}/${matrix}.mtx.gz) \
			-v | grep "time_mul" | cut -d' ' -f2 | tr -d '\n' >>gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt
			echo -n " " >>gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt
		done
	done
done

# end formats
for format in \
		"coo" \
		"csr" \
		"bsr -s 4" \
		"bsr -s 8" \
		"bsr -s 16"; do
		
		echo " " >>gp_$$_$(echo $format | tr ' ' '_').txt
		cat gp_$$_$(echo $format | tr ' ' '_').txt	>>gp_$$_res.txt
done
for i in 2 4; do
	for format in \
				"kat -s 4" \
				"kat -s 8" \
				"kat -s 16"; do
		ktmp="k=$i ${format}"
		echo " " >>gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt
		cat gp_$$_$(echo $ktmp | tr ' ' '_' | tr -d '"').txt >>gp_$$_res.txt
	done
done 

 >gp_$$_result.txt


exit 0
if false; then #----------------------------------------------------------------

echo -n "format "> gp_$$_header.txt
for matrix in $big_list; do echo "$matrix " | tr -d '\n' >gp_$$_${matrix}.txt; done

for format in \
		"coo" \
		"csr" \
		"bsr -s 16" \
		"bsr -s 32" \
		"bsr -s 64" \
		"bsr -s 128"; do
		
		echo -n "\"${format}\" " >>gp_$$_header.txt
done

time for matrix in $big_list; do

	for format in \
		"coo" \
		"csr" \
		"bsr -s 16" \
		"bsr -s 32" \
		"bsr -s 64" \
		"bsr -s 128"; do

		
		time ./main -f ${format} \
			-a <(gzip -cd ${big_dir_mat}/${matrix}.mtx.gz) \
			-b <(gzip -cd ${big_dir_vec}/vector_${matrix}_*.mtx.gz) \
			-V -v | grep "a_size" | cut -d' ' -f2 | tr -d '\n' >>gp_$$_${matrix}.txt
		
		echo -n " ">>gp_$$_${matrix}.txt
	done
done

for i in 2 4 8 16 32; do
	for format in \
				"kat -s 16" \
				"kat -s 32" \
				"kat -s 64" \
				"kat -s 128"; do
		echo -n "\"${format} k=$i\" " >>gp_$$_header.txt
	done
done 

for i in 2 4 8 16 32; do
	make DOUBLE_PRECISION=1 KAT_N=$i
	time for matrix in $big_list; do
		for format in \
			"kat -s 16" \
			"kat -s 32" \
			"kat -s 64" \
			"kat -s 128"; do
						
			time ./main -f ${format} \
				-a <(gzip -cd ${big_dir_mat}/${matrix}.mtx.gz) \
				-b <(gzip -cd ${big_dir_vec}/vector_${matrix}_*.mtx.gz) \
				-V -v -o ./resvec_${matrix}_$(echo ${format} | tr ' ' '_').mtx \
				| grep "a_size" | cut -d' ' -f2 | tr -d '\n' >>gp_$$_${matrix}.txt
			
			echo -n " ">>gp_$$_${matrix}.txt
		done
	done
done

echo " ">>gp_$$_header.txt
cat gp_$$_header.txt >gp_$$_result.txt

for matrix in $big_list; do
	echo " " >>gp_$$_${matrix}.txt
	cat gp_$$_${matrix}.txt >>gp_$$_result.txt
done

fi
