#!/bin/bash

#Description: This Shell Script iterates over various folders to read the
#information of various .lis fils. It must be placed in the directory containing
#all the folders of average strcuctures (1kf1_X-Yns_ave_struct). You need all
#.lis files ready within the folders.
#
# Argument #1: PDB ID (e.g. 1KF1, 1BNA)
#
#Requirements: Bash V4.0 or above
#Last Updated: 04-21-20 at 8:42pm

#Loop through all directories in the file

rm lislist
ls -v *.lis > lislist
LineCount=$(wc -l < lislist)
g++ -o WriteCSVH EfficientHReader.cpp
  ./WriteCSVH lislist $LineCount "$1.csv"
