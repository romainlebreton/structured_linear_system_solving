set xlabel "rank"
set ylabel "time"
set terminal pdf

set xrange [0:9499]
set out "crt.pdf"
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 1 lw 3 lc 4
plot "crtx10-before.dat" using 1:3 with linespoints title "before",\
     "crtx10-after.dat" using 1:3 with linespoints title "after",\
