#!/bin/bash

if [[ "$(hostname)" != "star" ]]; then
	source ./tests/utils.sh
else
	source /home/$(whoami)/eiaapp/tests/utils.sh
fi

LIMIT=1300

if [[ $# -eq 1 ]]; then
  LIMIT=$1
fi


echo "- loop_unrolling gcc_opt loop_unrolling+gcc_opt" > results

for i in $(seq 512 512 $LIMIT); do

	make_check "make"
A=$(valgrind --tool=cachegrind ./main mvm ./matrices/${i}.csr ./matrices/vector${i}.mtx 2>&1 \
		| grep "LL misses" | tr -s ' ' | cut -d ' ' -f 4 | tr -d ',' | tr '\n' ' ')

	make_check "MVLUR=1 make"
B=$(valgrind --tool=cachegrind ./main mvm ./matrices/${i}.csr ./matrices/vector${i}.mtx 2>&1 \
		| grep "LL misses" | tr -s ' ' | cut -d ' ' -f 4 | tr -d ',' | tr '\n' ' ')

	make_check "GCCOP=1 make"
C=$(valgrind --tool=cachegrind ./main mvm ./matrices/${i}.csr ./matrices/vector${i}.mtx 2>&1 \
		| grep "LL misses" | tr -s ' ' | cut -d ' ' -f 4 | tr -d ',' | tr '\n' ' ')

	make_check "GCCOP=1 MVLUR=1 make"
D=$(valgrind --tool=cachegrind ./main mvm ./matrices/${i}.csr ./matrices/vector${i}.mtx 2>&1 \
		| grep "LL misses" | tr -s ' ' | cut -d ' ' -f 4 | tr -d ',' | tr '\n' ' ')


MIN=$(echo $A $B $C $D | tr ' ' '\n' | sort -n | head -1)
A=$((A - MIN + 10))
B=$((B - MIN + 10))
C=$((C - MIN + 10))
D=$((D - MIN + 10))


	echo $A $B $C $D $i >> results

done

cat ./tests/cachegrind/gnuplot.conf | gnuplot

