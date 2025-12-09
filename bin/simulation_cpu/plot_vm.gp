# Set output terminal and filename for the plot
set terminal pngcairo size 800,600
set output 'results/Vm_Atrial_Control_VS_Drug.png'

# Set plot title and axis labels
set title 'Vm Grandi Model (Control VS Drug)'
set xlabel 'Time (ms)'
set ylabel 'Vm (millivolt)'
set datafile separator ","

# Set line style
set style line 1 linecolor rgb "black" linewidth 4
set style line 2 linecolor rgb "red" lt 0 linewidth 4
set style line 3 linecolor rgb "blue" lt 0 linewidth 4

# Plot data from the first file
# 'file1.dat' is the filename
# using 1:2 means use column 1 for x and column 2 for y
# with linespoints plots lines connecting points
# title 'Dataset 1' sets the legend entry for this plot
plot 'results/0.00/quinidine_0.00_time_series_smp0.csv' using 1:2 with lines linestyle 1 title 'Control', \
     'results/3237.00/quinidine_3237.00_time_series_smp0.csv' using 1:2 with lines linestyle 2 title 'Cmax1', \
     'results/6474.00/quinidine_6474.00_time_series_smp0.csv' using 1:2 with lines linestyle 3 title 'Cmax2'
