# PROCESS Long G4 STRUCTURE
#----------------------------------------------------------------------------
# This script works in conjunction with trajectory_smooth.tcl in order to get
# average structures of a trajectory file. This script  is specifically for one
# single trajectory file, $traj, defined by user. The file this was written for
# contains around 400,000 frames which is too big to load on to VMD. Thus, this
# script loads a window, $frames, defined by the user, and process that window for
# average structures.
#
# The Get_avg_pdbs procedure selects the Ab-Initio structures and saves them
# in a chosen directory. That procedure originates from the AVG_STRUCTS.tcl file.
#
# Arguments:
#   PDB - PDB ID for file names
#   PDBFile - This is the PDB structure to load for the trajectory
#   traj - The trajectory file
#
# Last Updated: 05/03/21
# By Angel-Emilio Villegas Sanchez

proc ProcessLongG4 {PDB PDBFile traj} {
    #Load the Averaging Script
    source trajectory_smooth.tcl

    #Create a directory for GG pairs to be stored
    file mkdir Pairs

    #Load PDB Structure
    mol new $PDBFile
    set iter 8
    set frames 50000
    set avgframes 1000

    #Create List of Trajectories to loop through
    set trajCount 0
    for {set i 0} {$i < $iter} {incr i} {
        animate delete all
        set FF [expr $i*$frames]
        set LF [expr $FF + [$frames - 1]]]
        mol addfile $traj first $FF last $LF waitfor all
        Get_avg_pdbs 0 $avgframes $PDB $trajCount
        set trajCount [expr $trajCount + 1]
    }
}

proc GGExtract {prefixname filenumber} {

set strands {{11 12 13 14} {34 35 36 37} {58 59 60 61}}

#Ths is the first index of the atoms we'll extract
set firstIndex 11

#This is the number of atoms we'll extract, minus 1
set atomNum 14

#This is the amount of pairs we'll extract with
set pairAmount 3

foreach residue $strands {
    for {set gg 0} {$gg < $pairAmount} {incr gg} {
        set res1 [lindex $residue $gg]
        set res2 [lindex $residue [expr $gg + 1]]
        set sel1 [atomselect top "residue $res1"]
        set sel2 [atomselect top "residue $res2"]

        #index that will be extracted
        set FirstG [lindex [$sel1 list] $firstIndex]
        set SecondG [lindex [$sel2 list] $firstIndex]
        $sel1 delete
        $sel2 delete

        #Make the selection and save to PDB
        set NextFG [expr $FirstG + $atomNum]
        set NextSG [expr $SecondG + $atomNum]

        #Make selection of 2 Guanines
        set sel [atomselect top "index $FirstG to $NextFG $SecondG to $NextSG"]

        #Create a string for naming
        set pairname g${res1}g${res2}

        $sel writepdb $prefixname\_$filenumber\_$pairname.pdb
        file copy $prefixname\_$filenumber\_$pairname.pdb Pairs/
        file delete $prefixname\_$filenumber\_$pairname.pdb
        $sel delete
    }
}

}

#In Directory GG Extraction
#--------------------------------------
#This procedure taking takes specific parameters to make
proc RangeGGExtract {prefixname step lim} {
    for {set i 0} {$i < $lim} {incr i} {
        set j [expr $i*100 + 49]
        mol top [mol new $prefixname\_$j.pdb]
        GGExtract $prefixname $j
        mol delete top
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
        set filenumber [expr 50*$trajCount + $i]
        avg_position $sel $prefixname\_$filenumber.pdb beg $BEG end $END writesel selonly

        # This conditional does an additional Ab-Initio Selection
        set AIFile [expr $i - 49]
        if {$AIFile == 0 || $AIFile % 100 == 0} {
            #Copy the file to another Folder
            file copy $prefixname\_$filenumber.pdb AbInitio/

            #Extract GGpair
            mol top [mol new $prefixname\_$filenumber.pdb]
            GGExtract $prefixname $filenumber
            mol delete top
        }
    }
}
