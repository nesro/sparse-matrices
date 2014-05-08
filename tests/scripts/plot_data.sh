#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

data_dir=${1:-../_out_mvm/out_mvm}

#-------------------------------------------------------------------------------
# warning: this script is ridiculously NOT efficient
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# outputs:
#  - tree_nodes.pdf
#  - mtx_sizes.pdf
#  - mul_speed.pdf
#-------------------------------------------------------------------------------

for i in $(ls $data_dir); do
	out=$(echo $i | cut -d '_' -f 1)
	method=$(echo $i | cut -d '_' -f 2)
	matrix_name=$(echo $i | cut -d '_' -f 3)
	matrix_format=$(echo $i | cut -d '_' -f 4)
	
	if [[ "$matrix_format" == "csr" ]]; then
		continue
	fi
	
	block_size=$(echo $i | cut -d '_' -f 5)
	kat_n=$(echo $i | cut -d '_' -f 6)
	mul_speed=$(grep "time_mul" $data_dir/$i | cut -d ' ' -f 2)
	#matrix_size=$(grep "a_n" $data_dir/$i | cut -d ' ' -f 2)
	matrix_size=666
	matrix_memory=$(grep "a_size" $data_dir/$i | cut -d ' ' -f 2)

	# We'll compare the speed with CSR
	csr_mul=$(ls ${data_dir}/*"$matrix_name"*"csr"* | xargs -i grep "time_mul" {} | cut -d ' ' -f 2)

	speedup_mul=$( bc <<< "scale=20; ${csr_mul}/${mul_speed}" | sed 's/^./0./' )
	speedup_ram=$( bc <<< "scale=20; (${matrix_size}*${matrix_size}*8)/${matrix_memory}" | sed 's/^./0./' )

	# if the matrix is in the kat format, save its tree properties
	if [[ "$matrix_format" == "kat" ]]; then
		kat_n_inner=$(grep "kat_a_inner" $data_dir/$i | cut -d ' ' -f 2)
		kat_n_dense=$(grep "kat_a_dense" $data_dir/$i | cut -d ' ' -f 2)
		kat_n_csr=$(grep "kat_a_csr" $data_dir/$i | cut -d ' ' -f 2)
		echo "$block_size-kat-$kat_n $kat_n_inner $kat_n_dense $kat_n_csr" \
			>>_$$_${matrix_name}_nodes.txt
	fi

	echo "out=$out method=$method name=$matrix_name format=$matrix_format bs=$block_size k=$kat_n \
mul_speed=$mul_speed size=$matrix_size ram=$matrix_memory speedupmul=$speedup_mul \
speedup_ram=$speedup_ram csr_mul=$csr_mul"
 
done

#-------------------------------------------------------------------------------

