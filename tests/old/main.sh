#!/bin/bash

source ./tests/utils.sh

SIZE=128

make

valgrind_check "./main generate ./tmp_matrices/vector$SIZE.mtx 1 $SIZE"
valgrind_check "./main generate ./tmp_matrices/$SIZE.mtx $SIZE $SIZE"
valgrind_check "./main convert ./tmp_matrices/$SIZE.mtx ./tmp_matrices/$SIZE.csr"
valgrind_check "./main mmm ./tmp_matrices/$SIZE.mtx"
valgrind_check "./main mmm ./tmp_matrices/$SIZE.csr"
valgrind_check "./main mvm ./tmp_matrices/$SIZE.mtx ./tmp_matrices/vector$SIZE.mtx"
valgrind_check "./main mvm ./tmp_matrices/$SIZE.csr ./tmp_matrices/vector$SIZE.mtx"


diff_check csr_mmm_hash.txt std_mmm_hash.txt
diff_check csr_mvm_hash.txt std_mvm_hash.txt

