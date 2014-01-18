#!/bin/bash

source ./tests/utils.sh

i=$(cat ./tests/start.nodel)

echo $((i + 512)) > ./tests/start.nodel

#echo "-O0+loop_unrolling -O3 -O3+loop_unrolling" > ./results.nodel

O0NO=
O0UR=
O3NO=
O3UR=

make_check "make clean && GCCOP=0 make"
TMP=$(./main eiammm $i foo)
O0NO=$(echo "$TMP" | grep -w "mmm_normal" | cut -d' ' -f2)
O0UR=$(echo "$TMP" | grep -w "mmm_unroll" | cut -d' ' -f2)
	
make_check "make clean && GCCOP=2 make"
TMP=$(./main eiammm $i foo)
O3NO=$(echo "$TMP" | grep -w "mmm_normal" | cut -d' ' -f2)
O3UR=$(echo "$TMP" | grep -w "mmm_unroll" | cut -d' ' -f2)

echo i=$i O0NO=$O0NO O0UR=$O0UR O3NO=$O3NO O3UR=$O3UR
echo i=$i O0NO=$O0NO O0UR=$O0UR O3NO=$O3NO O3UR=$O3UR >> ./tests/allresults.nodel
	
A="$(echo "$O0NO / $O0UR" | bc -l | sed 's/^\./0./')"
B="$(echo "$O0NO / $O3NO" | bc -l | sed 's/^\./0./')"
C="$(echo "$O0NO / $O3UR" | bc -l | sed 's/^\./0./')"

echo "$A $B $C $i"
echo "$A $B $C $i" >> ./tests/results.nodel
