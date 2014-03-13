#!/bin/bash

source ./tests/utils.sh

if [[ $# != 5 ]]; then
	echo "arguments: [mmm/mvm] limit_start limit_step limit_max threads"
	exit 1
fi

m=$1
limit_start=$2
limit_step=$3
limit_max=$4
threads=$5

filethreads=$(printf "%02d" $threads)

gpfile_time=./eia_${m}_gnuplot_time_${filethreads}.txt
gpfile_speedup=./eia_${m}_gnuplot_speedup_${filethreads}.txt
gpfile_png_time=eia_${m}_time_${filethreads}.png
gpfile_png_speedup=eia_${m}_speedup_${filethreads}.png

binary=./main

make_check "GCCOP=2 OMP_THREADS=$threads make"

sleep 1

methods="csr_${m} csr_${m}_unrolled csr_${m}_parallel csr_${m}_unrolled_parallel"

echo $methods >$gpfile_time 
echo $methods | cut -d' ' -f2- >$gpfile_speedup

for limit in $(seq $limit_start $limit_step $limit_max); do
	echo "$limit / $limit_max "
	for method in $methods; do
		echo -e "\tmethod=$method"
		eval "$method=\"$($binary eia $method $limit $(($limit * $limit / 2)))\""
		eval "echo -n \"\$$method \" >>$gpfile_time"
	done
	
	eval "echo -n \"\$(echo \" \$csr_${m} / \$csr_${m}_unrolled \" | bc -lq | sed 's/^\\./0./' ) \" >>$gpfile_speedup"
	eval "echo -n \"\$(echo \" \$csr_${m} / \$csr_${m}_parallel \" | bc -lq | sed 's/^\\./0./' ) \" >>$gpfile_speedup"
	eval "echo -n \"\$(echo \" \$csr_${m} / \$csr_${m}_unrolled_parallel \" | bc -lq | sed 's/^\\./0./' ) \" >>$gpfile_speedup"
			
	echo $limit >>$gpfile_time
	echo $limit >>$gpfile_speedup
done

gnuplot_wrapper $gpfile_time $gpfile_png_time "order of matrix" "time[s]" 5
gnuplot_wrapper $gpfile_speedup $gpfile_png_speedup "order of matrix" "speedup" 4
