#!/bin/bash

source ./tests/utils.sh

gpfile=./gnuplot_temp.txt

echo "first second third fourth" >$gpfile
echo "1 2 3 4 10" >>$gpfile
echo "2 3 4 5 20" >>$gpfile
echo "3 4 5 6 30" >>$gpfile

gnuplot_wrapper $gpfile ./gnuplot_temp.png xlabel ylabel 5