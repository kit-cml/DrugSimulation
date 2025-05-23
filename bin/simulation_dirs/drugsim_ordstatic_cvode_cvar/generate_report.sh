#!/bin/sh

if [ "$#" -lt 4 ]; then
    echo "Error: Provide result directory, number of samples, username, and LaTEX files sequentially"
    exit 1
fi

#Plot all the time-series result from the in-silico simulation
echo Generate plots from the time-series result
python3 ../../scripts/plot_time_series.py $1 $2

#Concat the separated feature data
echo Unifying feature data
python3 ../../scripts/plot_features.py $1 $3

#Generate report based on the pre-generated LaTEX file
echo "Generate PDF from LaTEX (on construction)"
pdflatex $4
