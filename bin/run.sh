#!/bin/bash

mass=$1
scale=`echo "$mass/2." | bc`
scaleDo=`echo "$mass/4." | bc`
scaleUp=$mass

dirce=${mass}_ce
dirup=${mass}_up
dirdo=${mass}_do

mkdir $dirce
mkdir $dirup
mkdir $dirdo

for i in `seq 1 100`; do 
  bsub -q cmscaf1nd -o /tmp/job_out "submit.sh $mass $scale   $(( 1000 + $i )) ${dirce}" 
  bsub -q cmscaf1nd -o /tmp/job_out "submit.sh $mass $scaleDo $(( 2000 + $i )) ${dirdo}" 
  bsub -q cmscaf1nd -o /tmp/job_out "submit.sh $mass $scaleUp $(( 3000 + $i )) ${dirup}"
done  
