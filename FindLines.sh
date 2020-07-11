#!/bin/bash
# This script is meant to go through every file and 
# check for a doubly existing line.

ls -v *.out > OutList

x=1
for thing in $(cat OutList)
do
	if [ 2 == $(grep -o "Thank you" $thing | wc -l) ]
	then
		echo "$x":Succesfully Complete!
	else 
		echo "$x":"$thing": Incomplete :{	
	fi
	((x++))	

done
