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

#Create CSV file
MyCSV="${1}_InterBase_Values.csv"
touch $MyCSV

#Loop through all directories in the file
for d in $(find $PWD -maxdepth 1 -type d)
do
  #Move CSV and Access directory
  mv $MyCSV $d
  cd $d

  #Execute Curves

  #Create LisList in directory to iterate through it
  ls -v *.lis > LisList

  for lis in 'cat $LisList' do

  done

  #Get out of folder
  cd ..
  pwd
  echo "Im ouutta here!"

done

#ls -v *.lis > LisList
