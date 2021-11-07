#Description: This script will help with the transformations of
proc align {auid ausel dnaid dnasel} {
     set dna [atomselect $dnaid "all"]
     set dnacensel [atomselect $dnaid $dnasel]
     set dnacenter [measure center $dnacensel]
     set dnaoffset [vecscale -1 $dnacenter]
     $dna moveby $dnaoffset
     set Au [atomselect $auid "all"]
     set Aucensel [atomselect $auid $ausel]
     set aucenter [measure center $Aucensel]
     set auoffset [vecscale -1 $aucenter]
     $Au moveby $auoffset

}


#Transform for AUII AUSS face
# $dna moveby $offset:
#Transform for AUI AUS Corner:
# $dna moveby $offset
# $dna move [transaxis y 45]
#Transform for AUI AUS Face:
# $dna moveby $offset <-- depends on distance to set the closest atom
# $dna move [transaxis y 54.702]
# $dna move [transaxis z 45]
#Transform for AUI AUS Face (45deg):
# $dna moveby $offset <-- offset to origin in order to rotate
# $dna move [transaxis x 45]
# $dna moveby $offset <-- to the distance outside the nanoparticle
# $dna move [transaxis y 54.702]
# $dna move [transaxis z 45]
