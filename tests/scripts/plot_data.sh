#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

#set -x

data_dir=${1:-../measured_data/create_graphs}
matrices=${2:-"EX6 fp gupta3 exdata-1 heart1  human-gene2"}

#gupta3

#-------------------------------------------------------------------------------
# warning: this script is ridiculously NOT efficient
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# outputs:
#  - tree_nodes.pdf
#  - mtx_sizes.pdf
#  - mul_speed.pdf
#-------------------------------------------------------------------------------

mkdir -p _plots_$$
tmpdir=$_

#-------------------------------------------------------------------------------

#prepare headers
echo -n "format " >$tmpdir/mmm-speedup.txt
echo -n "format " >$tmpdir/mvm-speedup.txt
echo -n "format " >$tmpdir/ram-up.txt

for matrix in $matrices; do
	echo -n "$matrix " >>$tmpdir/mmm-speedup.txt
	echo -n "$matrix " >>$tmpdir/mvm-speedup.txt
	echo -n "$matrix " >>$tmpdir/ram-up.txt
done

echo " " >>$tmpdir/mmm-speedup.txt
echo " " >>$tmpdir/mvm-speedup.txt
echo " " >>$tmpdir/ram-up.txt

#-------------------------------------------------------------------------------


for f in coo csr; do

	echo -n "\"$f\" " >>$tmpdir/ram-up.txt
	
	if [[ "$f" != "csr" ]]; then
	echo -n "\"$f\" " >>$tmpdir/mvm-speedup.txt
	fi
	
	for matrix in $matrices; do

#kat preparation
	echo "format vnitřní husté CSR" >$tmpdir/kat_nodes_${matrix}.txt

		echo "matrix $matrix"

		#ram
		den_size=$(  ls $data_dir/*mvm*$matrix*$f* | xargs -i grep a_n {} | awk '{ print $2 }' | tr -d '\n' )
		this_size=$( ls $data_dir/*mvm*$matrix*$f* | xargs -i grep a_size {} | awk '{ print $2 }' | tr -d '\n' )
		echo -n "$( bc <<< "scale=20; (${this_size}/(${den_size}*${den_size}*8))" ) " >>$tmpdir/ram-up.txt
		
		# Let's skip CSR
		if [[ "$f" == "csr" ]]; then
			continue;
		fi
		
		echo "mvm"
		
		#mvm
		csr_mvm_speed=$(  ls $data_dir/*mvm*$matrix*csr* | xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
		this_mvm_speed=$( ls $data_dir/*mvm*$matrix*$f*  | xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
		echo -n "$( bc <<< "scale=20; ${csr_mvm_speed}/${this_mvm_speed}" ) " >>$tmpdir/mvm-speedup.txt

	done
	
	echo " " >>$tmpdir/mvm-speedup.txt
	echo " " >>$tmpdir/ram-up.txt
done

	

#16 512
for bs in 32 64 128 256; do
	for f in bsr; do
	
		echo -n "\"${f}-${bs}\" " >>$tmpdir/ram-up.txt
		echo -n "\"${f}-${bs}\" " >>$tmpdir/mvm-speedup.txt
		echo -n "\"${f}-${bs}\" " >>$tmpdir/mmm-speedup.txt
	
		for matrix in $matrices; do
			#ram
			den_size=$(  ls $data_dir/*mvm*$matrix*${f}_${bs}* | xargs -i grep a_n {} | awk '{ print $2 }' | tr -d '\n' )
			this_size=$( ls $data_dir/*mvm*$matrix*${f}_${bs}* | xargs -i grep a_size {} | awk '{ print $2 }' | tr -d '\n' )
			echo -n "$( bc <<< "scale=20; (${this_size}/(${den_size}*${den_size}*8))" ) " >>$tmpdir/ram-up.txt
	
			#mvm
			csr_mvm_speed=$(  ls $data_dir/*mvm*$matrix*csr*        | xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
			this_mvm_speed=$( ls $data_dir/*mvm*$matrix*${f}_${bs}* | xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
			echo -n "$( bc <<< "scale=20; ${csr_mvm_speed}/${this_mvm_speed}" ) " >>$tmpdir/mvm-speedup.txt
			
			#mmm
			csr_mvm_speed=$(  ls $data_dir/*mmm*$matrix*csr*        | xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
			this_mvm_speed=$( ls $data_dir/*mmm*$matrix*${f}_${bs}* | xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
			
			if [[ "x$this_mvm_speed" == "x" ]]; then
				echo "emptty $f $bs";
			fi
			
			echo -n "$( bc <<< "scale=20; ${csr_mvm_speed}/${this_mvm_speed}" ) " >>$tmpdir/mmm-speedup.txt
	
		done
		
		echo " " >>$tmpdir/ram-up.txt
		echo " " >>$tmpdir/mvm-speedup.txt
		echo " " >>$tmpdir/mmm-speedup.txt
		
	done
	for f in kat; do
		
		for k in 2 4 8; do
		
			echo -n "\"${k}${f}-${bs}\" " >>$tmpdir/ram-up.txt
			echo -n "\"${k}${f}-${bs}\" " >>$tmpdir/mvm-speedup.txt
			echo -n "\"${k}${f}-${bs}\" " >>$tmpdir/mmm-speedup.txt
			echo -n "\"${k}${f}-${bs}\" " >>$tmpdir/kat-nodes.txt
	
			for matrix in $matrices; do
				#ram
				den_size=$(  ls $data_dir/*mvm*$matrix*${f}_${bs}_${k}* | xargs -i grep a_n {} | awk '{ print $2 }' | tr -d '\n' )
				this_size=$( ls $data_dir/*mvm*$matrix*${f}_${bs}_${k}* | xargs -i grep a_size {} | awk '{ print $2 }' | tr -d '\n' )
				echo -n "$( bc <<< "scale=20; (${this_size}/(${den_size}*${den_size}*8))" ) " >>$tmpdir/ram-up.txt
		
				#mvm
				csr_mvm_speed=$(  ls $data_dir/*mvm*$matrix*csr*             | \
					xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
				this_mvm_speed=$( ls $data_dir/*mvm*$matrix*${f}_${bs}_${k}* | \
					xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
				echo -n "$( bc <<< "scale=20; ${csr_mvm_speed}/${this_mvm_speed}" ) " >>$tmpdir/mvm-speedup.txt
			
				#mmm
				csr_mvm_speed=$(  ls $data_dir/*mmm*$matrix*csr*             | \
					xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
				this_mvm_speed=$( ls $data_dir/*mmm*$matrix*${f}_${bs}_${k}* | \
					xargs -i grep time_mul {} | awk '{ print $2 }' | tr -d '\n' )
					
				if [[ "x$this_mvm_speed" == "x" ]]; then
					echo "emptty $f $k $bs";
				fi
					
				echo -n "$( bc <<< "scale=20; ${csr_mvm_speed}/${this_mvm_speed}" ) " >>$tmpdir/mmm-speedup.txt
	
				#kat_nodes
	
	
	kat_n_inner=$( ls $data_dir/*mvm*$matrix*${f}_${bs}_${k}* | xargs -i grep kat_a_inner {} | awk '{ print $2 }' | tr -d '\n' )
	kat_n_dense=$( ls $data_dir/*mvm*$matrix*${f}_${bs}_${k}* | xargs -i grep kat_a_dense {} | awk '{ print $2 }' | tr -d '\n' )
	kat_n_csr=$( ls $data_dir/*mvm*$matrix*${f}_${bs}_${k}* | xargs -i grep kat_a_csr {} | awk '{ print $2 }' | tr -d '\n' )
		
				echo "\"${k}${f}-${bs}\" $kat_n_inner $kat_n_dense $kat_n_csr" \
					>>$tmpdir/kat_nodes_${matrix}.txt
	
			done
		
			echo " " >>$tmpdir/ram-up.txt
			echo " " >>$tmpdir/mvm-speedup.txt
			echo " " >>$tmpdir/mmm-speedup.txt
		done
	done
done





 ./tests/scripts/gnuplot.sh $tmpdir/mmm-speedup.txt 8 0 "formát matice" "zrychlení oproti CSR"
 ./tests/scripts/gnuplot.sh $tmpdir/mvm-speedup.txt 8 0 "formát matice" "zrychlení oproti CSR" 1
 ./tests/scripts/gnuplot.sh $tmpdir/ram-up.txt 8 1 "formát matice" "% velikost husté matice"
 
 ls  $tmpdir/kat_nodes* | xargs -i ./tests/scripts/gnuplot.sh {} 8 1 "formát matice" "počet uzlů"
 




exit 0


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

for bs in 16 32 64 128 256 512; do
	for f in bsr; do
		$script $matrix_dir $matrix $mul_vector $f $bs 2
		wait_for_slot
	done
	for f in kat; do
		for k in 2 4 8; do
			$script $matrix_dir $matrix $mul_vector $f $bs $k
			wait_for_slot
		done
	done
done


for mul_vector in 0 1; do
	for matrix in $matrices; do
		for f in coo csr; do
		
			# Let's skip SpMMM in the COO format
			if [[ "$f" == "coo" && "$mul_vector" == "0" ]]; then
				continue;
			fi
		
			$script $matrix_dir $matrix $mul_vector $f 0 2
			wait_for_slot
		done
		for bs in 16 32 64 128 256 512; do
			for f in bsr; do
				$script $matrix_dir $matrix $mul_vector $f $bs 2
				wait_for_slot
			done
			for f in kat; do
				for k in 2 4 8; do
					$script $matrix_dir $matrix $mul_vector $f $bs $k
					wait_for_slot
				done
			done
		done
	done
done

#-------------------------------------------------------------------------------

exit 0

#-------------------------------------------------------------------------------

echo -n "format " >$tmpdir/matrix_names.txt

for i in $(ls $data_dir); do
	out=$(echo $i | cut -d '_' -f 1)
	method=$(echo $i | cut -d '_' -f 2)
	matrix_name=$(echo $i | cut -d '_' -f 3)
	matrix_format=$(echo $i | cut -d '_' -f 4)
	
	block_size=$(echo $i | cut -d '_' -f 5)
	kat_n=$(echo $i | cut -d '_' -f 6)
	mul_speed=$(grep "time_mul" $data_dir/$i | cut -d ' ' -f 2)
	matrix_size=$(grep "a_n" $data_dir/$i | cut -d ' ' -f 2)
	matrix_size=666
	matrix_memory=$(grep "a_size" $data_dir/$i | cut -d ' ' -f 2)

	# We'll compare the speed with CSR
	csr_mul=$(ls ${data_dir}/*"${method}"*"$matrix_name"*"csr"* | xargs -i grep "time_mul" {} | cut -d ' ' -f 2)

	speedup_mul=$( bc <<< "scale=20; ${csr_mul}/${mul_speed}" )
	speedup_ram=$( bc <<< "scale=20; (${matrix_memory}/(${matrix_size}*${matrix_size}*8))" )

	# if the matrix is in the kat format, save its tree properties
	if [[ "$matrix_format" == "kat" && "$method" == "mmm" ]]; then
		kat_n_inner=$(grep "kat_a_inner" $data_dir/$i | cut -d ' ' -f 2)
		kat_n_dense=$(grep "kat_a_dense" $data_dir/$i | cut -d ' ' -f 2)
		kat_n_csr=$(grep "kat_a_csr" $data_dir/$i | cut -d ' ' -f 2)
		
		echo "\"kat.n=$kat_n bs=$block_size\" $kat_n_inner $kat_n_dense $kat_n_csr" \
			>>$tmpdir/_${matrix_name}_nodes.txt
	fi


	if ! grep $matrix_name $tmpdir/matrix_names.txt >/dev/null; then
		echo -n "$matrix_name" >>$tmpdir/matrix_names.txt
	fi

	if [[ "$matrix_format" == "kat" ]]; then
		echo "\"kat.n= $kat_n bs= $block_size $matrix_format\" $speedup_ram" >>$tmpdir/_${method}_${matrix_name}_ramup.txt
	elif [[ "$matrix_format" == "bsr" ]]; then
		echo "\"bs= $block_size $matrix_format\" $speedup_ram" >>$tmpdir/_${method}_${matrix_name}_ramup.txt
	else
		echo "\"$matrix_format\" $speedup_ram" >>$tmpdir/_${method}_${matrix_name}_ramup.txt	
	fi

	if [[ "$matrix_format" == "csr" ]]; then
		continue
	fi

	if [[ "$matrix_format" == "kat" ]]; then
		echo "\"kat.n= $kat_n bs= $block_size $matrix_format\" $speedup_mul" >>$tmpdir/_${method}_${matrix_name}_speedup.txt
	else
		echo "\"bs= $block_size $matrix_format\" $speedup_mul" >>$tmpdir/_${method}_${matrix_name}_speedup.txt	
	fi

done

cat $tmpdir/matrix_names.txt > $tmpdir/final.txt
ls $tmpdir | xargs -i sort $tmpdir/{} -n -k 2 -k 4 -o $tmpdir/{}
ls $tmpdir | xargs -i sed -i 's/= /=/g' $tmpdir/{}
ls $tmpdir | xargs -i sh -c "cat $tmpdir/matrix_names.txt >$tmpdir/final_{}; cat $tmpdir/{} >> $tmpdir/final_{}"


ls $tmpdir/final_* | xargs -i ./tests/scripts/gnuplot.sh {} 5

#-------------------------------------------------------------------------------

#	echo "out=$out method=$method name=$matrix_name format=$matrix_format bs=$block_size k=$kat_n \
#mul_speed=$mul_speed size=$matrix_size ram=$matrix_memory speedupmul=$speedup_mul \
#speedup_ram=$speedup_ram csr_mul=$csr_mul"s
