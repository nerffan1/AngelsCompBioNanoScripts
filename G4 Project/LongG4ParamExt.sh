#!/bin/bash

#Description: This is a quicker method to extract the parameters from the LIS files using Linux tools. The previous version was a custom-made C++ reader, but this script is less of a hassle/ 

#Get a list of files to process
if [ ! -f lislist ]
then
ls -v *.lis > lislist
fi

if [ ! -f LongG4Params ]
then
touch LongG4Params
fi

for lis in $(cat pairlist)
do
#Extract the areas within |G| and |H|
sed -n '/|G|/,/|H|/p' $lis | sed -n '/\(12)\)\|\(13)\)\|\(14)\)/'p | awk '{print $5,$6,$7,$8,$9,$10}' >> LongG4ParamsTI
done
