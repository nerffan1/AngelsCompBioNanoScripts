#!/bin/bash
#This program works for the submission of QCHEM input in the USC HPC.
#It's meant to Chiaka Onukwugha

#Function Description:
function CreateFromTags () {
for tag in $@
do 
   ls -v *$tag* >> TIList
done

sed -i 's/_input.inp//' TIList

num=0

for file in $(cat TIList)
do 
   cp base.slurm qc$num.slurm
   sed -i "s/inpo/$file/" qc$num.slurm
   sed -i "s/oup/$file.out/" qc$num.slurm
   sbatch qc$num.slurm
   (( num++ ))
done
rm TIList
}

#Function Description:
function Reconvergence () {
echo "Type in the file with the list of tags to reconverge:"
read Clean

echo Your current directory is $(pwd)   

echo "Where are the original input files located in relation to current directory?"
read OgLoc

#Find the location of the input files, in relation to your curent directory, and 

# cp ../GG_TI_Angel/*"$NewTag"* ./NEWINPUT        

num=1

for tag in $(cat $Clean)
do
	echo Copying ${OgLoc}${tag}_input.inp
	cp ${OgLoc}${tag}_input.inp .
	sed '43i scf_guess = read' ${tag}_input.inp > ${tag}.mid  
	sed -e '35r fix.txt' ${tag}.mid > ${tag}.inp

#The following section 
	cp base.slurm qc$num.slurm
   	sed -i "s/inpo/$tag.inp/" qc$num.slurm
   	sed -i "s/oup/$tag.out/" qc$num.slurm
	sbatch qc$num.slurm
   	(( num++ ))	
done 

}

function CreateCleanList(){
sed -e 's/.out//g' IncompleteList > NEWLISTCLEAN
        mkdir NEWINPUT
}

function RangeInput(){
num=0
for tag in $(eval echo {$1..$2..$3})
do
	ls -v *$tag* >> TILIST
done
sed -i 's/_input.inp//' TILIST

for file in $(cat TILIST)
do 
   cp base.slurm qc$num.slurm
   sed -i "s/inpo/$file/" qc$num.slurm
   sed -i "s/oup/$file.out/" qc$num.slurm
   sbatch qc$num.slurm
   (( num++ ))
done
rm -f TILIST
}

#Introduction:
cd /scratch2/angelemv/TESTDIR/
echo Hello, what would you like to do?
echo 1- Create QCHEM input from various tags
echo 2- Create QCHEM input for Reconvergence
echo 2a - Create Clean List
echo 3- Create individual QCHEM input
echo 4- Ranged QCHEM input
read choice1

case $choice1 in 
	1)
	echo -e "\nYou've chosen various tags to make your QCHEM input"	
	echo "\nPlease Enter your tags: "
	read input
	CreateFromTags ${input[@]}
	;;

	2)
	echo -e "\nYou've chosen reconvergence"
	Reconvergence
	;;
	
	2a)
	echo -e "\nCreate Clean list for Reconvergence from Output files"
	CreateCleanList
	;;

	3)
	echo -e  "\nYou've chosen to create individual QCHEM inputs"
	;;

	4) echo -e "\nSubmit QCHEM input giving initial, final, and increment: "
	read input
	RangeInput ${input[@]}	
	;;
esac

