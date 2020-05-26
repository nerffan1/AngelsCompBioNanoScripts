# Loop through every trajectory, and save them.
#
#Last Updated: 05/25/20

proc ProcessG4 {PDB} {
    set TrajList [glob *.trr]
    #source AVG_STRUCT.tcl
    mol new
    foreach TRR $TrajList {
        puts $TRR
    }
}
