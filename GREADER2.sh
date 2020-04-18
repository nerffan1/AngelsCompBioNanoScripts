#!/bin/bash

#Description: This Shell Script iterates over various folders to read the
#information of various .lis fils. It must be placed in the directory containing
#all the folders of average strcuctures (1kf1_X-Yns_ave_struct). You need all
#.lis files ready within the folders.
#
# Argument #1: PDB ID (e.g. 1KF1, 1BNA)
#
#Requirements: Bash V4.0 or above
#Last Updated: 04-03-20 at 8:13pm

#Loop through all directories in the file
for d in $(find . -maxdepth 1 -name "*_struct" -type d)
do
  #Access directory, create csv, and create crv/lis files
  cd $d
  touch MYCSV.csv
  ls -v *.pdb > pdblist
  source generateCRV.sh pdblist
  source curves.sh pdblist
  #Iterate through LisList
  ls -v *.lis > Lislist
  #List=$(find . -maxdepth 1 -name Lis* -type f)
  #echo $List
  #for lisfile in 'cat $List'
  #do
    source GREADERCSV $lisfile MYCSV.csv
  #done
  break
  cd ..
done

#ls -v *.lis > LisList
