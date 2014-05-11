#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

#-------------------------------------------------------------------------------
# parse the configuration.txt file
#-------------------------------------------------------------------------------


todo:



export username=nesrotom

#-------------------------------------------------------------------------------

export utils_debug=1

function debug_print {
	if [[ "$utils_debug" == "1" ]]; then
		echo "$@"
	fi
}

#-------------------------------------------------------------------------------

function gnuplot_bars {

	file_in=$1
	file_out=$2
	to=$2

	plot="plot '${file_in}.txt' using 2:xtic(1) ti col axes x1y2"
	for (( i=3; i<$to; i++ )); do
			plot+=", \"\" u $i ti col axes x1y2 "
	done  

	gnuplot <<__EOF__
set term pdf monochrome
set output '${file_out}.pdf'
set boxwidth 0.9 absolute
set style fill solid 1.00 border lt -1
set key inside right top vertical Left noreverse noenhanced autotitles nobox
#set style histogram clustered gap 1 title  offset character 0, 0, 0
set style fill pattern
set datafile missing '-'
set style data histograms
set xlabel "matrix"
set ylabel "time [s]"
set xtics border in scale 0,0 nomirror rotate by -45 offset character 0, 0, 0 autojustify
set xtics norangelimit font ",8"
set xtics ()
set logscale y2
unset ytics
set y2tics mirror
$plot
__EOF__

}

#-------------------------------------------------------------------------------

# We don't want to bother others. So we'll have 2 jobs max.
function wait_for_slot {

	if [[ "$(hostname)" != "star" ]]; then
		debug_print "You're not the STAR server!"
		return
	fi

	sleep 5
	while :; do
		if (( $( qstat | grep $username | wc -l ) < 2 )); then
			return
		fi
		sleep 5
	done
}

#-------------------------------------------------------------------------------

if false; then
	export big_list="cage12 ldoor Freescale1 thermal2 atmosmodj nlpkkt120"
	export big_dir_mat=../big_matrices
	export big_dir_vec=../big_vectors
else
	export big_list="c-23 eurqsa GT01R"
	export big_dir_mat=../test_matrices
	export big_dir_vec=../test_matrices
	
	export generated_dir=../test_matrices/generated
	export generated_list="test1 test2 test3 test4"
	
fi

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
	
	# set term png size 1024,786
	
	gnuplot << __EOF__
		set term pdf
		set output "$output_filename"
		set pointsize 2
		set style data point
		
		set logscale y

		set key top left reverse Left
		set key autotitle columnhead
		set key title "Legend"
		set key box width 1 height 1

		set ylabel "$ylabel"
		set xlabel "$xlabel"

		$plot
__EOF__
}

function check_code {
	code=$1
	if [[ $code -ne 0 ]]; then
		echo "Previous command has failed. Exiting."
		exit 1
	fi
}


#gnuplot_wrapper test.txt test.png "xlabel" "ylabel" 3
