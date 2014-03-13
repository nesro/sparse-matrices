#!/bin/bash

export RED='\e[0;31m'
export GREEN='\e[1;32m'
export YELLOW='\e[1;33m'
export PURPLE='\e[1;33m'
export LIGHT_PURPLE='\e[1;35m'
export NO_COLOR='\e[0m'

function make_check {
	echo -ne "make check: ${YELLOW}$1${NO_COLOR} "
	eval $1 >/dev/null
	if [[ "$?" == "0" ]] ; then
		echo -e "${GREEN}[OK]${NO_COLOR}"
	else
		echo -e "${RED}[FAIL]${NO_COLOR}"
		exit
	fi
}

function valgrind_check {
	echo -ne "valgrind check: ${YELLOW}$1${NO_COLOR} "
	valgrind --leak-check=full --show-reachable=yes $1 2>&1 | grep "ERROR SUMMARY: 0 errors from 0 contexts" >/dev/null
	if [[ "$?" == "0" ]] ; then
		echo -e "${GREEN}[OK]${NO_COLOR}"
	else
		echo -e "${RED}[FAIL]${NO_COLOR}"
	fi
}

function diff_check {
	echo -n "diff check: $1 $2 "
	if [[ "$(diff $1 $2)" == "" ]] ; then
		echo -e "${GREEN}[OK]${NO_COLOR}"
	else
		echo -e "${RED}[FAIL]${NO_COLOR}"
	fi
}

function run_command {
	echo -e "${LIGHT_PURPLE}[running:${YELLOW} $1${LIGHT_PURPLE}]${NO_COLOR} "
	echo " V "
	eval $1
	echo " ^ "
}

function gnuplot_wrapper {
	input_filename=$1
	output_filename=$2
	xlabel=$3
	ylabel=$4
	nlines=$5
	
	plot="plot \"$input_filename\" using $nlines:1"
	for (( i=2; i<nlines; i++ )); do
		plot+=", \"\" using $nlines:$i"
	done  
	
	gnuplot << __EOF__
		set term png size 1024,786
		set output "$output_filename"
		set style data linespoints

		set key top left reverse Left
		set key autotitle columnhead
		set key title "Legend"
		set key box width 1 height 1

		set ylabel "$ylabel"
		set xlabel "$xlabel"

		$plot
__EOF__
}

