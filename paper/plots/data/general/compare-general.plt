set xlabel "size"
set ylabel "time in sec."
set terminal pdf
set size ratio 2

set out "compare-general.pdf"
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 1 lw 8 lc 8
set key left top
set xtics 1100

set multiplot 
set title "5 blocks with n columns"
set size ratio 1.3
set yrange [0:4500]
set origin -0.45,0
plot "testx5-dac.dat" using 1:3 with lines ls 1 title "dac",\
     "testx5-newton.dat" using 1:3 with lines ls 3 title "newton",\
     "testx5-crt.dat" using 1:3 with lines ls 2 title "crt"
set title "n blocks with n columns"
set xtics 3000
set ylabel ""
set size ratio 1.3
set origin 0.25,0
set yrange [0:50000]
set xrange [0:9000]
plot "dac.dat" using 1:3 with lines ls 1 title "dac",\
     "crt.dat" using 1:3 with lines ls 2 title "crt"

unset multiplot
