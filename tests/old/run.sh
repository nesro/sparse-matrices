#!/bin/bash

TEST_DIR=./tests
TEST_BIN=./tests/bin
VALGRIND_RESULTS=./tests/valgrind

function run_test {
	test_binary=$1
	valgrind_log=${VALGRIND_RESULTS}/valgrind_$(date +%d.%m.%Y_%H:%M:%S).txt
	
	valgrind --log-file=$valgrind_log --leak-check=full --show-reachable=yes \
		$test_binary

	if ! grep -q "ERROR SUMMARY: 0 errors from 0 contexts" $valgrind_log; then
		echo "valgrind error"
		cat $valgrind_log
	fi
}

for i in $(ls $TEST_BIN); do
	run_test $TEST_BIN/$i
done
