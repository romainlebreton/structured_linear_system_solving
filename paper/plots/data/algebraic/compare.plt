set xlabel "size"
set ylabel "time"
set terminal pdf
set size ratio 2

set out "compare-algebraic.pdf"
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 0 lw 8 lc 8

set multiplot 
set title "5 blocks with n columns"
set size ratio 1.3
set xrange [0:10000]
set origin -0.45,0
plot "testx5-bmc.dat" using 1:3 with linespoints title "bmc",\
     "testx5-naive.dat" using 1:3 with linespoints title "naive"

set title "n blocks with n columns"
set ylabel ""
set size ratio 1.3
set xrange [0:10000]
set origin 0.25,0
plot "test1000-bmc.dat" using 1:3 with linespoints title "bmc",\
     "test1000-naive.dat" using 1:3 with lines ls 2 title "naive"

unset multiplot
