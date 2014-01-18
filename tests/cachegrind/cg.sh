#!/bin/bash

if [[ $# -eq 1 ]]; then
  LIMIT=$1
fi

echo "-O0+loop_unrolling -O3 -O3+loop_unrolling" > ./tests/cachegrind/results.txt

START=512
STEP=512
LIMIT=600

for i in $(seq $START $STEP $LIMIT); do
	echo i=$i

	make clean 2>&1 >/dev/null
	DEBUG=1 make 2>&1 >/dev/null
	valgrind --tool=cachegrind ./main cn $i a 2>&1 >/dev/null
	#Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
	totals=$(callgrind_annotate --auto=yes cachegrind* | grep "PROGRAM TOTALS" | tr -d ',')
	O0NO=$(echo $totals | cut -d' ' -f 6)
	
	make clean 2>&1 >/dev/null
	GCCOP=2 DEBUG=1 make 2>&1 >/dev/null
	valgrind --tool=cachegrind ./main cn $i a 2>&1 >/dev/null
	#Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
	totals=$(callgrind_annotate --auto=yes cachegrind* | grep "PROGRAM TOTALS" | tr -d ',')
	O0UR=$(echo $totals | cut -d' ' -f 6)
	
	make clean 2>&1 >/dev/null
	DEBUG=1 make 2>&1 >/dev/null
	valgrind --tool=cachegrind ./main cu $i a 2>&1 >/dev/null
	#Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
	totals=$(callgrind_annotate --auto=yes cachegrind* | grep "PROGRAM TOTALS" | tr -d ',')
	O3NO=$(echo $totals | cut -d' ' -f 6)
	
	make clean 2>&1 >/dev/null
	GCCOP=2 DEBUG=1 make 2>&1 >/dev/null
	valgrind --tool=cachegrind ./main cu $i a 2>&1 >/dev/null
	#Ir I1mr ILmr Dr D1mr DLmr Dw D1mw DLmw
	totals=$(callgrind_annotate --auto=yes cachegrind* | grep "PROGRAM TOTALS" | tr -d ',')
	O3UR=$(echo $totals | cut -d' ' -f 6)


	A="$(echo "$O0NO / $O0UR" | bc -l | sed 's/^\./0./')"
	B="$(echo "$O0NO / $O3NO" | bc -l | sed 's/^\./0./')"
	C="$(echo "$O0NO / $O3UR" | bc -l | sed 's/^\./0./')"

	echo "$A $B $C $i"
	
	echo "$A $B $C $i" >> ./tests/cachegrind/results.txt
	
done

(echo "set term png size 1024,786" && echo "set output \"cache_grind.png\"" && \
	cat ./tests/cachegrind/gnuplot.conf) | gnuplot

