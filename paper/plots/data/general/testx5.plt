set xlabel "rank"
set ylabel "time"
set terminal pdf

set xrange [0:5000]
set yrange [0:4500]
set out "testx5-general.pdf"
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 1 lw 3 lc 4
plot "testx5-dac.dat" using 1:3 with linespoints title "dac",\
     "testx5-newton.dat" using 1:3 with linespoints title "newton",\
     "testx5-dixon.dat" using 1:3 with linespoints title "dixon",\
     "testx5-crt.dat" using 1:3 with linespoints title "crt",\

