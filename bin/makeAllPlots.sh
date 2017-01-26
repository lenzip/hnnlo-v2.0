#!/bin/bash

for i in `ls | grep "_ce" | grep -v PtJetMin50`; do
  basename=`echo $i | cut -d "_" -f 1`
  python makePlot.py $basename
done
