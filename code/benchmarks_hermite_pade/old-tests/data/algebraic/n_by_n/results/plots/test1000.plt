set xlabel "rank"
set ylabel "time"
set terminal pdf

set out "test1000-algebraic.pdf"
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 1 lw 3 lc 4
plot "test1000-bmc.dat" using 1:3 with linespoints title "bmc",\
     "test1000-naive.dat" using 1:3 with linespoints ls 2 title "naive"
