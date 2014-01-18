#!/bin/bash

if [[ "$(hostname)" != "star" ]]; then
	source ./tests/utils.sh
else
	source /home/$(whoami)/eiaapp/tests/utils.sh
fi

# ---------------------------------------------------------------------------

#echo "MMM TEST"

#MATRIX_SIZE=1024

#rm -f "csr_mmm_hash.txt" "unroll_csr_mmm_hash.txt"

#make_check "GCCOP=2 make"

#run_command "./main mmm ./matrices/$MATRIX_SIZE.csr"

#make_check "GCCOP=2 MMLUR=1 make"

#run_command "./main mmm ./matrices/$MATRIX_SIZE.csr"

#diff_check "csr_mmm_hash.txt" "unroll_csr_mmm_hash.txt"

# ---------------------------------------------------------------------------

echo "MVM TEST"

if [[ $# -gt 2 ]]; then
  if [ "$3" -eq "$3" ] 2> /dev/null; then
    MATRIX_SIZE=$3
  else
    echo "not a number"
    exit
  fi
else
  MATRIX_SIZE=4096
fi

#make clean
make

MATRIX_MTX="matrices/$MATRIX_SIZE.mtx"
MATRIX_CSR="matrices/$MATRIX_SIZE.csr"
VECTOR="matrices/vector$MATRIX_SIZE.mtx"

touch tmp_matrices/necum
if [ ! -f matrices/$MATRIX_SIZE.mtx ]; then
  echo "file matrices/$MATRIX_SIZE.mtx not found!"
  echo "generate tmp_matrices/$MATRIX_SIZE.mtx ..."
  ./main generate tmp_matrices/$MATRIX_SIZE.mtx $MATRIX_SIZE $MATRIX_SIZE
  MATRIX_MTX="tmp_matrices/$MATRIX_SIZE.mtx"
  echo "done!"
fi

if [ ! -f matrices/$MATRIX_SIZE.csr ]; then
  echo "file matrices/$MATRIX_SIZE.csr not found!"
  echo "convert $MATRIX_MTX to tmp_matrices/$MATRIX_SIZE.csr ..."
  echo "./main convert $MATRIX_MTX tmp_matrices/$MATRIX_SIZE.csr"
  ./main convert $MATRIX_MTX tmp_matrices/$MATRIX_SIZE.csr
  MATRIX_CSR="tmp_matrices/$MATRIX_SIZE.csr"
  echo done!
fi

if [ ! -f matrices/vector$MATRIX_SIZE.mtx ]; then
  echo "file matrices/vector$MATRIX_SIZE.mtx not found!"
  echo "generate tmp_matrices/vector$MATRIX_SIZE.mtx ..."
  ./main generate tmp_matrices/vector$MATRIX_SIZE.mtx 1 $MATRIX_SIZE
  VECTOR="tmp_matrices/vector$MATRIX_SIZE.mtx"
  echo "done!"
fi

rm -f "csr_mvm_hash.txt" "unroll_csr_mvm_hash.txt"

rm *o
unset MVLUR
unset GCCOP
make_check "make"

run_command "./main mvm $MATRIX_CSR $VECTOR"

rm *o
make_check "MVLUR=2 make"

run_command "./main mvm $MATRIX_CSR $VECTOR"

rm *o
unset MVLUR
make_check "GCCOP=2 make"

run_command "./main mvm $MATRIX_CSR $VECTOR"

rm *o
make_check "GCCOP=2 MVLUR=1 make"

run_command "./main mvm $MATRIX_CSR $VECTOR"

diff_check "csr_mvm_hash.txt" "unroll_csr_mvm_hash.txt"

# ---------------------------------------------------------------------------
