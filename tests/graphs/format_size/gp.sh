#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

gp=./tests/graphs/format_size/gp

plot="plot '$gp.txt' using 2:xtic(1) ti col axes x1y2"
for (( i=3; i<6; i++ )); do
		plot+=", \"\" u $i ti col axes x1y2"
done  

gnuplot <<__EOF__
	set term pdf monochrome
	set output '$gp.pdf'
	set boxwidth 0.9 absolute
	set style fill solid 1.00 border lt -1
	set key inside right top vertical Left noreverse noenhanced autotitles nobox
	#set style histogram clustered gap 1 title  offset character 0, 0, 0
	set datafile missing '-'
	set style data histograms
	set style fill pattern
	set xlabel "matrix"
	set ylabel "size [B]"
	set xtics border in scale 0,0 nomirror rotate by -45 offset character 0, 0, 0 autojustify
	set xtics norangelimit font ",8"
	set xtics ()
	set logscale y2
	unset ytics
	set y2tics mirror
	$plot
__EOF__

#w: plot '$gp.txt' using 2:xtic(1) ti col axes x1y2, '' u 3 ti col axes x1y2, "" u 4 ti col axes x1y2, "" u 5 ti col axes x1y2

# plot '$gp.txt' using 2:xtic(1) ti col axes x1y2, '' u 3 ti col axes x1y2, \
	#using 4:xtic(1) ti col axes x1y2, '' u 4 ti col axes x1y2