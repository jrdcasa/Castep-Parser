reset
set style line 1 lt 1 ps 1.0 lc rgb "black"  pt 6 lw 1.0
set style line 2 lt 1 ps 1.0 lc rgb "red"    pt 5 lw 1.0
set style line 3 lt 2 ps 0.4 lc rgb "blue"   pt 4 lw 2.0
set style line 4 lt 1 ps 0.4 lc rgb "green"  pt 4 lw 2.0
set style line 5 lt 2 ps 0.4 lc rgb "yellow" pt 4 lw 2.0
set style line 6 lt 2 ps 0.4 lc rgb "orange" pt 4 lw 2.0

###############################################################################
set term wxt 1 enhanced dashed size 1200,800 font "Arial,10"
set multiplot layout 4,4

filelist = system("ls ./02-ITER_DAT/*.dat")
set title noenhanced
do for [file in filelist] {
  set title sprintf("%s", file)
  p file u 1:2 w p ls 1 notitle, file u 1:2 w l ls 1 notitle
}

unset multiplot

