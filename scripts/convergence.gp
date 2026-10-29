# run from the repository root: gnuplot scripts/convergence.gp
set terminal pngcairo size 1000,500
set output 'convergence/convergence.png'
set multiplot layout 1,2

set title 'L2 self-convergence error vs time'
set xlabel 't'
set ylabel '|| Pi_n - Pi_2n ||'
set logscale y
set key bottom right
plot 'convergence/error.dat' index 0 using 1:3 with linespoints title 'dr = pi/500 (n=250 vs 500)', \
     'convergence/error.dat' index 1 using 1:3 with linespoints title 'dr = pi/1000 (n=500 vs 1000)', \
     'convergence/error.dat' index 2 using 1:3 with linespoints title 'dr = pi/2000 (n=1000 vs 2000)'

unset logscale y
set title 'Measured convergence order'
set ylabel 'order'
set yrange [0:8]
plot 'convergence/order.dat' using 1:2 with linespoints title '250/500/1000', \
     'convergence/order.dat' using 1:3 with linespoints title '500/1000/2000', \
     4 with lines dt 2 title '4th order'
unset multiplot
