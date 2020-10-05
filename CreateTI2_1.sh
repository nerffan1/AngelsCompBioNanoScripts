#!/bin/bash
#This program works for the submission of QCHEM input in the USC HPC.
#It has several functions to submit files to the HPC

#Function Description:
function CreateFromTags () {
cd Input/
for tag in $@
do
	ls -v *$tag* >> TILIST
done
mv TILIST ..
cd ..
sed -i 's/_input.inp//' TILIST
num=0
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

#Function Description:
function Reconvergence () {
echo "Make sure to perform this function in the directory holding Input and Output directories!"
echo "Type in the file with the list of clean tags to reconverge:"
read Clean

echo Your current directory is $(pwd)   

echo "Where are the original input files located in relation to current directory? (include '/') "
read OgLoc

#Find the location of the input files, in relation to your curent directory, and 

# cp ../GG_TI_Angel/*"$NewTag"* ./NEWINPUT        

num=1

#Makes a directory called ReconvergeFiles 
echo "A directory called 'ReconvergeFiles' will be created if not found"
mkdir -p ReconvergeFiles
cp fix.txt ReconvergeFiles
echo "The input files to be reconverged will be copied from your Input directory to Reconverge Files"
#Ask user if he wishes to coninue
echo "Do you wish to continue?"
read cont
if [[ "$cont" != "y" ]]
then
	exit 1
fi 
#Loop through clean input list. Copy input files 
for tag in $(cat $Clean)
do
	echo Copying ${OgLoc}${tag}_input.inp to ReconvergeFiles
	cp ${OgLoc}${tag}_input.inp ReconvergeFiles
	cd ReconvergeFiles
	sed '43i scf_guess = read' ${tag}_input.inp > ${tag}.mid  
	sed -e '35r fix.txt' ${tag}.mid > ${tag}.inp
	rm ${tag}.mid
#The following section 
	cp base.slurm qc$num.slurm
   	SetUpBatch $tag $num
   	(( num++ ))	
done 

}

#This should be ran in the output section.
function CreateCleanList(){
sed -e 's/.out//g' IncompleteList > NEWLISTCLEAN
}
#Description: Rename base.slurm
#Parameters: (1) Tag, (2) Slurm Batch Number
function SetUpBatch(){
cp base.slurm qc$2.slurm
sed -i "s/inpo/$1/" qc$2.slurm
sed -i "s/oup/$1.out/" qc$2.slurm
sbatch qc$2.slurm
 
}

function RangeInput(){
num=0
cd Input/
for tag in $(eval echo {$1..$2..$3})
do
	ls -v *$tag* >> TILIST
done
mv TILIST ..
cd ..
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
echo 5- Resubmit Unsuccesful Output
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

	5) echo -e "\nSubmit unsuccesful Output List "
	read file
	sed -i 's/.out//' $file
	CreateFromTags $(cat $file)
	;;
esac

