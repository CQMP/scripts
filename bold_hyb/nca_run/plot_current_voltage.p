set terminal pdf enhanced font ",18" size 9,4

set output "output/nca_U2.0ed-0.1bw10.pdf"

set key left top
set yrange[0:*]
set xlabel "V, meV" offset 0,0.5
set ylabel "I" offset 1,0

set for [i = 1:10] style line i ps 1.5 lc i+1
set style line 5 ps 1.2 lc rgb "orange"

temp_vals="0.07 0.1 0.15 0.2 0.3 0.5 1.0"
T(n) = word(temp_vals,n)

set multiplot
set origin 0,0
set size 0.5,1
plot \
    for [i = 1:7] 'nca_stats/current_U2.0ed-0.1bw10.0T'.T(i).'t5.0.dat' u 1:2:3 w yerrorbars t 'T='.T(i) ls i

set origin 0.5,0
set yrange [*:*]
set ylabel "dI/dV" offset 2,0
set xlabel "T, meV" offset 0,0.5
set nokey
plot \
    'nca_stats/cond_t2.0.dat' u 1:2 w lp lc 2 lw 2 pt 6

