#!/bin/bash

if (( $# == 2 )); then
	n=$1
	s=$2
	d=$3
else
	n=8
	s=2
	d=50
fi

let max_items=$(($n*$n*($d/100)))
let items=0

echo "%%MatrixMarket matrix coordinate real general"
echo "$n $n $(($n * $n))"

for (( i=1; i <= $n; i++ )); do
	for (( j=1; j <= $n; j++ )); do
		echo $i $j $s
		(( s += 1 ))
	done
done
