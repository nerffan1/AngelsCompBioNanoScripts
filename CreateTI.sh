#!/bin/bash

for tag in $@
do
   ls -v *$tag* >> TIList
done

sed -i 's/_input.inp//' TIList

num=1

for file in $(cat TIList)
do
   cp qc.slurm qc$num.slurm
   sed -i "s/inp/$file.inp/" qc$num.slurm
   sed -i "s/oup/$file.out/" qc$num.slurm
   sbatch qc$num.slurm
   (( num++ ))
done
rm TIList
