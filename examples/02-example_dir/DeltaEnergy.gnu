reset
set style line 1 lt 1 ps 1.0 lc rgb "red"    pt 5 lw 1.0
set style line 2 lt 1 ps 1.0 lc rgb "black"  pt 6 lw 1.0
set style line 3 lt 2 ps 0.4 lc rgb "blue"   pt 4 lw 2.0
set style line 4 lt 1 ps 0.4 lc rgb "green"  pt 4 lw 2.0
set style line 5 lt 2 ps 0.4 lc rgb "yellow" pt 4 lw 2.0
set style line 6 lt 2 ps 0.4 lc rgb "orange" pt 4 lw 2.0

###############################################################################
set term wxt 1 enhanced dashed size 400,400 font "Arial,8"
set multiplot layout 1,1
set encoding utf8

f1="./02-ITER_DAT/Energy.dat"

set ylabel "{/Symbol D}E_{rel} (kcal/mol)"
set format y  "%.1f"
set key left top
set bmargin at screen 0.25


set xtics 1.0
set mxtics 2
set mytics 2
set xtics rotate by 90 right
set xtics noenhanced

p f1 u ($5==1?$4:$4/0):xtic(1) w p ls 1 title "Not Converged", f1 u 4:xtic(1) w l ls 1 notitle,\
  f1 u ($5==0?$4:$4/0):xtic(1) w p ls 2 title "Converged", f1 u 4:xtic(1) w l ls 2 notitle

unset multiplot

