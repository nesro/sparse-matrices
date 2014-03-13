#!/bin/bash

DIR=./tests/eia/mvm

source $DIR/../../utils.sh

START=8092
STEP=1024
LIMIT=15000

if [[ $# -eq 1 ]]; then
  LIMIT=$1
fi

echo "-O0+loop_unrolling -O3 -O3+loop_unrolling" > ./tests/results.txt

echo START=$START
echo STEP=$STEP
echo LIMIT=$LIMIT

O0NO=
O0UR=
O3NO=
O3UR=

for i in $(seq $START $STEP $LIMIT); do
	echo i=$i
	
	# FIXME: staci 2x na celeej
	
	make_check "make clean && GCCOP=0 make"
	TMP=$(./main eiamvm $i foo)
	O0NO=$(echo "$TMP" | grep -w "mvm_normal" | cut -d' ' -f2)
	O0UR=$(echo "$TMP" | grep -w "mvm_unroll" | cut -d' ' -f2)
	
	make_check "make clean && GCCOP=2 make"
	TMP=$(./main eiamvm $i foo)
	O3NO=$(echo "$TMP" | grep -w "mvm_normal" | cut -d' ' -f2)
	O3UR=$(echo "$TMP" | grep -w "mvm_unroll" | cut -d' ' -f2)
	
	A="$(echo "$O0NO / $O0UR" | bc -l | sed 's/^\./0./')"
	B="$(echo "$O0NO / $O3NO" | bc -l | sed 's/^\./0./')"
	C="$(echo "$O0NO / $O3UR" | bc -l | sed 's/^\./0./')"
	
	echo "$A $B $C $i" >> ./tests/results.txt
	
done

(echo "set term png size 1024,786" && echo "set output \"mvm.png\"" && \
	cat ./tests/gnuplot.conf) | gnuplot
