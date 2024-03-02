set terminal png size 640,480
set output 'strong_plot.png'
set title 'time to complete (1 core)'
set xrange [:1025]
set xlabel 'cores'
set ylabel 'sec'
set xtics (64,128,256,512,1024)
set logscale x 2
set logscale y 2
set key right bottom
plot for [col=2:2] 'strong.dat' using 1:col with lines title columnheader lw 2
