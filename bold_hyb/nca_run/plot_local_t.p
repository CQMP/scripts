set terminal pdf enhanced font ",18" size 9,4

set output "output/n_t_U2.0ed-0.1bw10.pdf"

set key left top
set xlabel "t, 1/meV" offset 0,0.5
set ylabel "n" offset 1,0
set yrange [0.5:1.0]

set for [i = 1:10] style line i ps 1.5 lc i+1
set style line 5 ps 1.2 lc rgb "orange"

temp_vals="0.07 0.1 0.15 0.2 0.3 0.5 1.0"
T(n) = word(temp_vals,n)

set multiplot
set origin 0,0
set size 0.5,1

set key maxcols 2 maxrows 4

V="0.05"

set label 1 "V=".V." meV" font ",20" right at graph 0.95, 0.05

plot \
    for [i = 1:7] 'data/G0.3/T'.T(i).'/bw10.0/U2.0/ed-0.1/v'.V.'/tmax5.0/n0_t.dat' u 1:($2*2) w l t 'T='.T(i) ls i

set origin 0.5, 0
V="0.2" 
set label 1 "V=".V." meV" font ",20" right at graph 0.95, 0.05
plot \
    for [i = 1:7] 'data/G0.3/T'.T(i).'/bw10.0/U2.0/ed-0.1/v'.V.'/tmax5.0/n0_t.dat' u 1:($2*2) w l t 'T='.T(i) ls i

