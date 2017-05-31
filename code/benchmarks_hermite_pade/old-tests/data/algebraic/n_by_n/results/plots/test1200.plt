set xlabel "rank"
set ylabel "time"
set terminal pdf

set out "test1200-algebraic.pdf"
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 1 lw 3 lc 4
plot "test1200-bmc.dat" using 1:3 with linespoints title "bmc",\
     "test1200-naive.dat" using 1:3 with lines ls 2 title "naive"
