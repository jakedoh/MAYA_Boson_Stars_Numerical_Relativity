set size 0.6
set nologscale y
set title '8-4 operators'
set xlabel 'Time'
set ylabel 'error (L_{/Symbol \245})'
set key right bottom

p "data/sevenpatches-plane-8-4-err-hi3/error.maximum.asc" u 2:3 title 'dx=0.0125 (min. bandwidth)' w l lw 3
rep "data/sevenpatches-plane-8-4-err-min-err-hi3/error.maximum.asc" u 2:3 title 'dx=0.0125 (min. ABTE)' w l lw 3

set out 'comparison-8-4-color.eps'
set terminal postscript eps enhanced color
replot
set terminal x11
