#!/bin/bash
for i in {0..1500000..5000}
#for i in {0..400..100}
do
	echo $i
	echo "29" | gmx_mpi trjconv -s ../../Cas9_2_eq5_RNAprotDNA.gro  -f  ../../traj_every50ps_DNA-prot-RNA_cas9_corr.xtc -o ${i}_ps.pdb -dump ${i} -n ../index.ndx

       sed -i 's/ DA/  A/g;s/ DG/  G/g;s/ DC/  C/g;s/ DT/  T/g' ${i}_ps.pdb
            # changing the bases
    
       cp ${i}_ps.pdb common.pdb
         Cur5 < curves_inp
       mv  common.lis ${i}_ps.lis
       mv  common_grp.pdb ${i}_ps_grp.pdb
done


