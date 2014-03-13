#!/bin/bash
set -x
for (( i=256; i<=16384; i*=2 )); do
	for (( j=2; j<=3; j++ )); do
		time ./gen ~/mm_mtx/${i}_0$((j-1)).mtx $i $j  
	done
done

