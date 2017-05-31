set xlabel "size"
set ylabel "time"
set terminal pdf
set size ratio 2

set out "compare-general.pdf"
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 0 lw 8 lc 8
set xtics 0,1000,5000

set multiplot 
set title "n by small"
set size ratio 1.3
set yrange [0:4500]
set origin -0.45,0
plot "testx5-dac.dat" using 1:3 with linespoints title "dac",\
     "testx5-newton.dat" using 1:3 with linespoints title "newton",\
     "testx5-dixon.dat" using 1:3 with linespoints title "dixon",\
     "testx5-crt.dat" using 1:3 with linespoints title "crt",\

set title "n by n"
set ylabel ""
set size ratio 1.3
set origin 0.25,0
set xrange [0:5500]
set yrange [0:50000]
plot "crt.dat" using 1:3 with linespoints title "crt",\
     "dac.dat" using 1:3 with linespoints title "dac",\

unset multiplot
