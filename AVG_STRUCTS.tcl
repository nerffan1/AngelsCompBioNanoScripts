#AVERAGE STRUCTURES OVER A TRAJECTORY
#--------------------------------------
#This code seeks to find the average structure at
#various intervals in a trajectory.
#MAKE SURE TO Source, or load, the trajectory_smooth.tcl script

proc Get_avg_pdbs {molid nstructs prefixname} {
    set TotalFrames [molinfo $molid get numframes]
    if {$TotalFrames % $nstructs != 0} {
        puts "Your frames divided by your structures must have 0 remainder."
        return 0
    }
    set FramesPS [expr $TotalFrames/$nstructs]
    set sel [atomselect $molid "nucleic"]

    #Loop through the frames and create files for each
    for {set i 0} {$i < $nstructs} {incr i} {
        set BEG [expr $i*$FramesPS]
        set END [expr $BEG + $FramesPS - 1]
        avg_position $sel $prefixname\_$i.pdb beg $BEG end $END writesel selonly
    }
}
