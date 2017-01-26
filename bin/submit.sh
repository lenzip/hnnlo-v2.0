#!/bin/bash

home=/afs/cern.ch/user/l/lenzip/work/ww2016/HighMassMoriond2017/jetbinning/hnnlo-v2.0/bin/
card=HNNLO-LHC13.input

mass=$1
scale=$2
seed=$3
dir=$4

cp $home/hnnlo .
cp $home/PDFsets . -r
cp $home/Pdfdata . -r
cp $home/br* . 

cat ${home}/${card} | sed -e "s#MASS#$mass#g" | sed -e "s#SCALE#$scale#g" | sed -e "s#SEED#$seed#g" > input.input

cat input.input

./hnnlo < input.input | tee ${seed}.log

cp *.top $home/$dir
cp ${seed}.log $home/$dir




