#AVERAGE STRUCTURES OVER A TRAJECTORY
#--------------------------------------
#This code seeks to find the average structure at
#various intervals in a trajectory.
#MAKE SURE TO Source, or load, the trajectory_smooth.tcl script, along this script, in VMD.
# Arguments:
#   molid - The VMD Molecule ID
#   windown - The window number, or the number of frames to average by
#   prefixname - A prefix to attach to the output pdb files (e.g. 1kf1_9501-10000ns_ave_struct)
#Last Update: 04/09/20 at 09:54pm

proc Get_avg_pdbs {molid windown prefixname} {
    source trajectory_smooth.tcl
    set TotalFrames [molinfo $molid get numframes]
    if {$TotalFrames % $windown != 0} {
        puts "Your frames divided by your structures has no 0 remainder."
        puts "There'll be frames that are not averaged"
    }
    set nstructs [expr $TotalFrames/$windown]
    puts "The number of structures is $nstructs"
    set sel [atomselect $molid "nucleic"]

    #Loop through the frames and create files for each
    for {set i 0} {$i < $nstructs} {incr i} {
        set BEG [expr $i*$windown]
        set END [expr $BEG + $windown - 1]
        avg_position $sel $prefixname\_$i.pdb beg $BEG end $END writesel selonly
    }
}
