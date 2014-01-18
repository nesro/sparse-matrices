#!/bin/bash

if (( $# == 2 )); then
	n=$1
	s=$2
else
	n=8
	s=2
fi

echo "%%MatrixMarket matrix coordinate real general"
echo "$n $n $(($n * $n))"

for (( i=1; i <= $n; i++ )); do
	for (( j=1; j <= $n; j++ )); do
		echo $i $j $s
		(( s += 1 ))
	done
done
