#!/bin/bash
# Once you generate the crv file (see generateCRVforPDB.sh), you can actually run Curves program on pdbs
if [ $# -lt 1 ]
then
	echo -e "\n\tUsage: $0 <pdbList>"
	echo -e "\t\t<pdbList>: list of pdbs in a file for which you have to generate the lis files"
	echo -e "\n\texiting ..."
	exit
fi

pdbList=$1

for pdb in `cat $pdbList`
do

	name=`basename $pdb .pdb`
	crvfile="$name".crv
	sed -i back 's:DA: A:;s:DC: C:;s:DG: G:;s:DT: T:g' $pdb  && ~/Documents/structural_biology_software/Curves5.3/Cur5 < $crvfile && rm $pdb"back"
done
