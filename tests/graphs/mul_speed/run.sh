#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

set -x

source ./tests/scripts/utils.sh


retrieve=mul_speed

type=matvec
#type=matmat

bsron=0
bsrtyp="csr bsr_-s_16 bsr_-s_32 bsr_-s_64 bsr_-s_128 bsr_-s_256"

katon=1
kattyp="kat_-s_16 kat_-s_32 kat_-s_64 kat_-s_128 kat_-s_256"

generated_list=
generated_dir=
generated_dir_vectors=

#-------------------------------------------------------------------------------

echo -n "format "> gp_$$_res.txt
for matrix in $generated_list; do echo "$matrix " | tr -d '\n' >>gp_$$_res.txt; done
echo " ">>gp_$$_res.txt

# begin formats
if (( $bsron == 1 )); then
	for format in $bsrtyp; do
		formatTR=$(echo $format | tr '_' ' ')
		echo -n "\"${formatTR}\" " >gp_$$_$format.txt
	done
fi
if (( $katon == 1 )); then
	for i in 2 4; do
		for format in $kattyp; do
			formatTR=$(echo $format | tr '_' ' ')
			ktmp="k=$i ${formatTR}"
			echo -n "\"$ktmp\" " >>gp_$$_${i}_${format}.txt
		done
	done 
fi


if (( $bsron == 1 )); then
make PRECISION=2
for matrix in $generated_list; do
	for format in $bsrtyp; do
		formatTR=$(echo $format | tr '_' ' ')

		time ./main -f ${formatTR} \
			-a <(gzip -cd ${generated_dir}/${matrix}.mtx.gz) \
			-v | grep "$retrieve" | cut -d' ' -f2 | \
			tr -d '\n' >>gp_$$_$format.txt
		
		echo -n " ">>gp_$$_$format.txt
	done
done
fi

if (( $katon == 1 )); then
	for i in 2 4; do
		make PRECISION=2 KAT_N=$i
		time for matrix in $generated_list; do
			for format in $kattyp; do
				formatTR=$(echo $format | tr '_' ' ')
				ktmp="k=$i ${format}"
						
				time ./main -f ${formatTR} \
				-a <(gzip -cd ${generated_dir}/${matrix}.mtx.gz) \
				-v | grep "$retrieve" | cut -d' ' -f2 | \
				tr -d '\n' >>gp_$$_${i}_$format.txt
				echo -n " " >>gp_$$_${i}_$format.txt
			done
		done
	done
fi

# end formats
if (( $bsron == 1 )); then
	for format in $bsrtyp; do
		formatTR=$(echo $format | tr '_' ' ')
		echo " " >>gp_$$_$format.txt
		cat gp_$$_$format.txt >>gp_$$_res.txt
done
fi
if (( $katon == 1 )); then
	for i in 2 4; do
		for format in $kattyp; do
			formatTR=$(echo $format | tr '_' ' ')
			ktmp="k=$i ${formatTR}"
			echo " " >>gp_$$_${i}_${format}.txt
			cat gp_$$_${i}_${format}.txt >>gp_$$_res.txt
		done
	done 
fi

