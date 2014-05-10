#!/bin/bash
# Tomas Nesrovnal, nesro@nesro.cz, Copyright 2014
# https://github.com/nesro/sparse-matrices

file=${1%%.*}
to=$2

logscale=${3:-0}

xlabel=${4:-"xlabel"}
ylabel=${5:-"ylabel"}

keyoutside=${6:-0}

if [[ $logscale = 0 ]]; then
	logscale_str=""
else
	logscale_str="set logscale y2"
fi

if [[ $keyoutside = 0 ]]; then
	key_str="set key inside right top vertical Left noreverse noenhanced autotitles nobox"
else
	key_str="set key out vert center top Left noreverse noenhanced autotitles nobox"
fi



plot="plot '${file}.txt' using 2:xtic(1) ti col axes x1y2"
for (( i=3; i<$to; i++ )); do
		plot+=", \"\" u $i ti col axes x1y2 "
done  

	gnuplot <<__EOF__
set term pdf monochrome
set output '${file}.pdf'
set boxwidth 0.9 absolute
set style fill solid 1.00 border lt -1
#set key outside
#set key inside right top vertical Left noreverse noenhanced autotitles nobox
#set key out vert center top Left noreverse noenhanced autotitles nobox
#set style histogram clustered gap 1 title  offset character 0, 0, 0
$key_str
set style fill pattern
set datafile missing '-'
set style data histograms
set xlabel "$xlabel"
set ylabel "$ylabel"
set xtics border in scale 0,0 nomirror rotate by -45 offset character 0, 0, 0 autojustify
set xtics norangelimit font ",8"
set xtics ()
$logscale_str
unset ytics
set y2tics mirror
$plot
__EOF__
