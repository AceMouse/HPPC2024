set terminal png size 640,480
set output 'weak_plot.png'
set title 'time to complete (1024*1024 nfreq)'
set xrange [:129]
set xlabel 'cores'
set ylabel 'sec'
set xtics (1,2,4,16,32)
set logscale x 2
set logscale y 2
set key right bottom
plot for [col=2:2] 'weak.dat' using 1:col with lines title columnheader lw 2
