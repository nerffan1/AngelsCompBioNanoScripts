# PROCESS G4 STRUCTURES
#----------------------------------------------------------------------------
# This script works in conjunction with AVG_STRUCTS and
# with trajectory_smooth.tcl in order to do 3 tasks:
#   1) Get Average structures of a trajectory file
#   2) Get a selection of these structures for Ab-Initio Transfer-Integral calculations
#   3) Get a selection of Guanine pairs for each of those selected structures
# All trajectories are assumed to be in the .trr format in this script, but could
# work with other compatible formats as long as "mol addfile" accepts it.
#
# The Get_avg_pdbs procedure selects the Ab-Initio structures and saves them
# in a chosen directory. That procedure originates from the AVG_STRUCTS.tcl file.
#
# Arguments:
#   PDB - PDB ID for file names
#   PDBFile - This is the PDB structure to load for the trajectory
#
#Last Updated: 05/26/20
# By Angel-Emilio Villegas Sanchez

proc ProcessG4 {PDB PDBFile} {
    #Load the Averaging Script
    source trajectory_smooth.tcl

    #Load PDB Structure
    mol new $PDBFile

    #Create List of Trajectories to loop through
    set trajCount 0
    foreach traj [glob *.trr] {
        animate delete all
        mol addfile $traj waitfor all
        Get_avg_pdbs 0 100 $PDB $trajCount
        set trajCount [expr $trajCount + 1]
    }
}

#AVERAGE STRUCTURES OVER A TRAJECTORY + G Pair Extraction
#--------------------------------------
# This script finds the average structure at various intervals in a trajectory, but has
# an additional functionality of extracting Guanine Pairs. The Guanine extraction needs
# to be set up manually, by giving it the range of Guanine atoms, and the name of the pair
# that will be used in the filename.
#
# Arguments:
# molid - ID of molecule that has loaded trajectory
# windown - number of frames to average by
# prefixname - Name to prefix to the files
# trajCount - This helps the numbering of files
#

proc Get_avg_pdbs {molid windown prefixname trajCount} {
    set TotalFrames [molinfo $molid get numframes]
    if {$TotalFrames % $windown != 0} {
        puts "Your frames divided by your structures has no 0 remainder."
        puts "There'll be frames that are not averaged"
    }
    set nstructs [expr $TotalFrames/$windown]
    set sel [atomselect $molid "nucleic"]

    #Loop through the frames and create files for each
    for {set i 0} {$i < $nstructs} {incr i} {
        set BEG [expr $i*$windown]
        set END [expr $BEG + $windown - 1]
        set filenumber [expr 500*$trajCount + $i]
        avg_position $sel $prefixname\_$filenumber.pdb beg $BEG end $END writesel selonly

        # This conditional does an additional Ab-Initio Selection
        set AIFile [expr $i - 49]
        if {$AIFile == 0 || $AIFile % 100 == 0} {
            #Copy the file to another Folder
            puts "I will copy the Ab_initio file"
            file copy $prefixname\_$filenumber.pdb AbInitio/

            #Guanine Extraction
            foreach {f j k} {"74 to 88" "107 to 121" "g3g4" "107 to 121" "140 to 154" "g4g5"} {
                set sel2 [atomselect top "index $f $j"]
                $sel2 writepdb $prefixname\_$filenumber\_$k.pdb
                file copy $prefixname\_$filenumber\_$k.pdb Pairs/
                file delete $prefixname\_$filenumber\_$k.pdb
                $sel2 delete
            }
            foreach {f j k} {"269 to 283" "302 to 316" "g9g10" "302 to 316" "335 to 349" "g10g11"} {
                set sel2 [atomselect top "index $f $j"]
                $sel2 writepdb $prefixname\_$filenumber\_$k.pdb
                file copy $prefixname\_$filenumber\_$k.pdb Pairs/
                file delete $prefixname\_$filenumber\_$k.pdb
                $sel2 delete
            }
            foreach {f j k} {"463 to 477" "496 to 510" "g15g16" "496 to 510" "529 to 543" "g16g17"} {
                set sel2 [atomselect top "index $f $j"]
                $sel2 writepdb $prefixname\_$filenumber\_$k.pdb
                file copy $prefixname\_$filenumber\_$k.pdb Pairs/
                file delete $prefixname\_$filenumber\_$k.pdb
                $sel2 delete
            }
            foreach {f j k} {"658 to 672" "691 to 705" "g21g22" "691 to 705" "724 to 738" "g22g23"} {
                set sel2 [atomselect top "index $f $j"]
                $sel2 writepdb $prefixname\_$filenumber\_$k.pdb
                file copy $prefixname\_$filenumber\_$k.pdb Pairs/
                file delete $prefixname\_$filenumber\_$k.pdb
                $sel2 delete
            }
        }
    }
}
