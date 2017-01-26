#!/bin/bash
for i in `find -name total.top`; do
  rm $i; 
done  

for i in `ls | grep "_ce"`; do
  basename=`echo $i | cut -d "_" -f 1`
  cd ${basename}_ce
  ../merge.sh
  cd ../${basename}_up
  ../merge.sh
  cd ../${basename}_do
  ../merge.sh
  cd ..
done

