#!/bin/bash
final=""
echo $1/*.csv
for f in $1/*.csv; 
do 
echo $(basename $f)
./pcorrelation -u $f -o $(basename $f).pdf -t Correlation Plot; 
done
