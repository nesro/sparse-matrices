#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

if [[ $# -ne 1 ]]; then
	echo "usage: $0 binary]"
	exit 1
fi

if [[ ! -x $1 ]]; then
	echo "$0: $1 is not executable, exiting"
	exit 1
fi

nice -n 19 valgrind --leak-check=full --show-reachable=yes \
	$1 2> >(while read l; do echo -e "\e[01;31m$l\e[0m" >&2; done)
