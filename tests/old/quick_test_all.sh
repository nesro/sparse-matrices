#!/bin/bash

source ./tests/utils.sh
gpfile=./gnuplot_temp.txt
mtxdir=~/mm_mtx

lo=1024
hi=8192

if [[ $# == 1 && $1 == "gnuplot" ]] ; then
	gnuplot_wrapper $gpfile ./gnuplot_temp.png matrix_size time 5
	exit 0
fi

if [[ $# == 3 && $1 == "lohi" ]]; then
	lo=$2
	hi=$3
else
	echo "csr quadtree-128 quadtree-256 quadtree-512 quadtree-1024" >$gpfile
fi

set -x

for (( i=$lo; i<=$hi; i*=2 )); do
	res=$(./main -f csr -a $mtxdir/${i}_01.mtx -b $mtxdir/${i}_02.mtx -o none)
	echo "$res" | grep mul | cut -d' ' -f2 | tr -d '\n' >>$gpfile
	echo -n " " >>$gpfile

	for (( j=128; j<=1024; j*=2 )); do
		res=$(./main -f quadtree -l $j -a $mtxdir/${i}_01.mtx -b $mtxdir/${i}_02.mtx -o none)
		echo "$res" | grep mul | cut -d' ' -f2 | tr -d '\n' >>$gpfile
		echo -n " " >>$gpfile
	done

	echo $i >>$gpfile

done

gnuplot_wrapper $gpfile ./gnuplot_temp.png "matrix size" "time[s]" 6
