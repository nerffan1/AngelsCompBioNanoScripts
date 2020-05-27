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
# Arguments:
#   PDB - PDB ID for file names
#   PDBFile - This is the PDB structure to load for the trajectory
#Last Updated: 05/26/20
# By Angel-Emilio Villegas Sanchez

proc ProcessG4 {PDB PDBFile} {
    #Load the Averaging Script
    source trajectory_smooth.tcl

    #Load PDB Structure and delete frame 0
    mol new $PDBFile
    animate delete beg 0 end 0 skip 0 0

    #Create List of Trajectories to loop through
    foreach traj [glob *.trr] {
        mol addfile $traj
        Get_avg_pdbs 0 100 $PDB
    }

}

#AVERAGE STRUCTURES OVER A TRAJECTORY
#--------------------------------------
# This code seeks to find the average structure atvarious intervals in a trajectory.

proc Get_avg_pdbs {molid windown prefixname} {
    source trajectory_smooth.tcl
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
        avg_position $sel $prefixname\_$i.pdb beg $BEG end $END writesel selonly

        # This conditional does an additional Ab-Initio Selection
        set AIFile [expr $i - 49]
        if {$AIFile == 0 || $AIFile % 100 == 0} {
            # Copy the file to another Folder

            #Do the Guanine Pair extraction
        }
    }
}
