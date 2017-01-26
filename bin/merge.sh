#!/bin/bash

../mergedata 1 `ls *.top | tr "\n" " "`
mv fort.12 total.top
python ../gnuplot2root.py total.top
