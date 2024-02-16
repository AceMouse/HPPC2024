set terminal png size 640,480
set output 'perf_plot.png'
set title 'time to complete'
set xrange [:129]
set xlabel 'no mol'
set ylabel 'sec'
set xtics (2,4,16,128)
set logscale x 2
set logscale y 2
plot for [col=2:3] 'perf.dat' using 1:col with lines title columnheader

