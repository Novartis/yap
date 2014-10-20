#!/bin/bash
gnuplot << EOF
set terminal postscript eps enhanced color solid rounded 
set output "$2"
set xlabel "Read Position"
set ylabel "Base Count"
set title " Position vs. Base Count"
 plot "$1" using 1:2 with lines t "A", \
"" using 1:3 with lines t "C", \
""  using 1:4 with lines t "T", \
""  using 1:5 with lines t "G", \
""  using 1:6 with lines t "N"
EOF
