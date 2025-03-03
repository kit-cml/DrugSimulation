#!/usr/bin/env gnuplot
set terminal pngcairo font "Arial,14" size 1024,768
set output "vm_ord_cvode.png"
set title "Plot of Vm of EP Simulation of ORd in CVode"
set xlabel "Time (msec)"
set ylabel "Vm (mV)"
set datafile separator ","

plot 'result/time_series_marcell.plt' using 1:2 with lines lw 4 \
