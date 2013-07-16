#!/bin/bash
barcodef=barcodes
plots="aligned_split"
minimum=0
maximum=1000000
echo    "MAKE SURE TO CONFIGURE THIS SCRIPT BEFORE RUNNING IT"
echo    "current config: barcode:$barcodef plots:$plots"
for line in $(cat $barcodef); 
do 
    cd      $plots
    mkdir   $line
    mv      *$line* $line
    echo    "mv *$line* $line"
    cd      ..
done
./plotsc | bash
