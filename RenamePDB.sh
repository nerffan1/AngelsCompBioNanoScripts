#!/bin/bash
#Description: This script changes the previous naming convention of
# 1kf1_RANGEns_ave_struct_N.pdb to simply 1kf1
#Give first number as A parameters
ls -v *.pdb > LIST
(( ia = $1 ))
(( ib = $1 ))
for d in $(cat LIST)
do
echo $d
echo $ia
#mv $d "1kf1_${i}.pdb"
((ia++))
done
#Check Output
echo "Is this the proper Output?"
read conf
if [[ "$conf" == "y" ]];
then
  echo "Your change of filenames will begin shortly..."
  sleep 7
  for d in $(cat LIST)
  do
  mv $d "1kf1_${ib}.pdb"
  ((ib++))
  done
else
  echo "Your change of filenames have been aborted."
fi
rm LIST
#ls -v *.lis > LisList
