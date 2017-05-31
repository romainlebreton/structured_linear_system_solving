set xlabel "rank"
set ylabel "time"
set terminal pdf

set out "compare-nbits.pdf"
set style line 1 lt 1 lw 8 lc 4
set style line 2 lt 1 lw 4 lc 7
set style line 3 lt 1 lw 3 lc 4
plot "test500-bmc.dat" using 1:3 with linespoints title "nbits=500",\
     "test800-bmc.dat" using 1:3 with linespoints title "nbits=800",\
     "test900-bmc.dat" using 1:3 with linespoints title "nbits=900",\
     "test1000-bmc.dat" using 1:3 with linespoints title "nbits=1000",\
     "test1100-bmc.dat" using 1:3 with linespoints title "nbits=1100",\
     "test1200-bmc.dat" using 1:3 with linespoints title "nbits=1200",\
     "test1300-bmc.dat" using 1:3 with linespoints title "nbits=1300",\
     "test1400-bmc.dat" using 1:3 with linespoints title "nbits=1400",\
     "test1500-bmc.dat" using 1:3 with linespoints title "nbits=1500",\
     "test1800-bmc.dat" using 1:3 with linespoints title "nbits=1800",\
     "test1900-bmc.dat" using 1:3 with linespoints title "nbits=1900",\
     "test2000-bmc.dat" using 1:3 with linespoints title "nbits=2000",\
     
