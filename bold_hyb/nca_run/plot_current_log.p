set terminal pdf enhanced font ",18" size 9,4
V="0.05"
bw="5.0"
ed="-0.1"

set output "output/current_log_v".V.".pdf"

set key left top
#set yrange[0:*]
set xlabel "t, 1/meV" offset 0,0.5
set ylabel "(I - I(tmax)) / I(tmax)" offset 2,0

set for [i = 1:10] style line i ps 1.5 lc i+1
set style line 5 ps 1.2 lc rgb "orange"

temp_vals="0.03 0.05 0.07 0.1 0.15 0.2 0.3 0.5 1.0"
T(n) = word(temp_vals,n)
dataf(T,V)="nca_stats/current_t_U2.0ed".ed."bw".bw."T".T."v".V.".dat"

set multiplot
set origin 0,0
set size 0.64,1

# plot log-lin current
set key maxcols 2 maxrows 5 right top width -2 samplen 2
set yrange [0.0001:1.0]
set logscale y

plot \
    for [i = 1:7] dataf(T(i),V) u 1:2 w l t 'T='.T(i) ls i, \
    for [i = 1:7] dataf(T(i),V) u 1:3 w l notitle ls i dashtype 2 lw 2, \
     dataf(T(1),V) u 1:(1e2) w l t "exp (-{/Symbol a}(T) t)" ls 1 lc -1 dashtype 2 lw 2

# plot exp fit
set origin 0.6,0
set size 0.4,1
unset logscale
set yrange [*:*]
datafit(V)='nca_stats/cur_fit_U2.0ed'.ed.'bw'.bw.'v'.V.'.dat'
set xtics 0.2
set mxtics 5
set ylabel "{/Symbol a}"
set xlabel "T, meV"
plot \
    datafit(V) u 1:4:5 w yerrorbars ps 1.4 pt 6 lc rgb "navy-blue" t "{/Symbol a}", \
    '' u 1:($1/0.5)+1.0 w l t "2*T + 1"
