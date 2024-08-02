set terminal png
set output "gnuplot.representation.png"
set title "Amplicon Representation"
set style line 1 lt 1 lc rgb "orange"
set style line 2 lt 1 lc rgb "green"
set style fill solid
set style data histograms
set xrange [0: mx]

my=system("echo $ampstatsgnuplotmy")
set yrange [-my:my]

fdata=system("echo $ampstatsgnuplotfd")
rdata=system("echo $ampstatsgnuplotrd")
print fdata

plot fdata notitle, rdata notitle
