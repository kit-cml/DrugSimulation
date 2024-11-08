#!/bin/sh

#Plot all the time-series result from the in-silico simulation
echo Generate plots from the time-series result
python3 plotting.py

#Concat the separated feature data
echo Unifying feature data
python3 concat.py

#Generate report based on the pre-generated LaTEX file
echo Generate PDF from LaTEX
pdflatex report_drug.tex
