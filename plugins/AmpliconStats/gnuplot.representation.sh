#! /bin/sh 
maxx=$1
export ampstatsgnuplotmy=$2
export ampstatsgnuplotfd=$3
export ampstatsgnuplotrd=$4


/usr/bin/gnuplot -e "mx=$maxx" $5'gnuplot.representation.gnuplot'
