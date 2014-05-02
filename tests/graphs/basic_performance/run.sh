#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

set -x

source ./tests/scripts/utils.sh

prog=./main
generator=./tests/bin/matrix_generator
mtx_dir=./tests/matrices/generated
mm_dir=./tests/matrices/matrix_market
gp=./tests/graphs/basic_performance/gp
tmp_dir=/var/tmp/basic_performance
mat=${tmp_dir}/matrix.mtx

function run_type {
	format=$1
	info=$2 # time_mul|size
	file=$3
	
	if [[ $# -ge 4 ]]; then
		block_size="-s $4"
	else
		block_size=""
	fi

	res=$( $prog -a $mat -f $format -v $block_size )
	echo -n "$format " >>$file
	echo $res | grep "$info" | cut -d ' ' -f 2 | tr '\n' ' ' >>$file
	echo $res | grep "$info" | cut -d ' ' -f 2 | tr '\n' ' ' >>$file
	echo "" >>$file

}


if false; then
	make tests
	check_code $?
fi

if [[ ! -f "$mat" ]]; then
	mkdir -p $tmp_dir
	nice -n 19 $generator -n 1024 -m -s 0 -irblocks,2,0,64,64,0.95 -S 2 -o $mat
	check_code $?
fi

if false; then
	echo "format a b" >$gp.txt
	run_type coo time_mul $gp.txt
fi

gnuplot <<__EOF__
	set term pdf
	set output '$gp.pdf'
	set boxwidth 0.9 absolute
	set style fill solid 1.00 border lt -1
	set key inside right top vertical Left noreverse noenhanced autotitles nobox
	set style histogram clustered gap 1 title  offset character 0, 0, 0
	set datafile missing '-'
	set style data histograms
	set xlabel "corpus file"
	set ylabel "memory usage [MB]"
	set xtics border in scale 0,0 nomirror rotate by -45 offset character 0, 0, 0 autojustify
	set xtics norangelimit font ",8"
	set xtics ()
	#set logscale y2
	unset ytics
	set y2tics mirror
	plot '$gp.txt' using 2:xtic(1) ti col axes x1y2, '' u 3 ti col axes x1y2
__EOF__
