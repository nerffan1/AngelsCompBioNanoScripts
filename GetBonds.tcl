# GET GetGGBonds
#------------------------------------------------------------------
# Description: The GetGGBonds procedure gets the bonds from a set of files that
# are (1) located in the directory, and (2) have the same structure (i.e. the same
# amount of atoms and residues). The pair indices are obtained by processing
# the first file in an ordered list, $listfile, provided by the user. These
# pairs are printed in order for the user to save them.
# A set of bond lengths will be appended to a csv file. A

package require csv

proc GetGGBonds {listfile} {
     # Read a list since TCL orders files differently with globbing
     set rfile [open $listfile r]
     set files [read $rfile]
     close $rfile

     set first 1
     foreach struct $files {
          mol load pdb $struct
          set sel [atomselect top "all"]
          if {$first == 1} {
               set Pairs [GetBondPairs $sel]
               puts [list $Pairs] ;# Pairs are printed for user to save
               set first 0
          }
          WriteCSV $Pairs
          mol delete all
     }
}
#Description: This function takes in an atomselect selection in order to
proc GetBondPairs {sel} {
     set BondPairs {}
     set bonds [$sel getbonds]
     set siz [llength $bonds]
     #For every atom i
     for {set atomNum 0} {$atomNum < $siz} {incr atomNum} {
          set atom [lindex [$sel list] $atomNum]
          foreach bondedToSet [lindex $bonds $atomNum] { ;#For every set of atoms, s, bonded to i
               foreach bondedTo $bondedToSet { ;#For every atom in s
                    set recBond [lsort [list $atom $bondedTo]]
                    set check [lsearch $BondPairs $recBond]
                    if { $check == -1} { ;#Only append if it's not there
                         lappend BondPairs $recBond
                    }
               }
          }
     }
     return lsort -integer -index 0 $BondPairs ;#Order based on first value
}

# Description: After obtaining a list with all pair of indices corresponding to bonds, WriteCSV
# appends, to a CSV file, the bond values. This output is customized for a paticular PDB structure
# containing 2 Guanine structures.
proc WriteCSV {Pairs} {
set file [open bondlengths.csv a+]
for {set pair 0} {$pair < 16} {incr pair} {
     lappend bondl [measure bond [lindex $Pairs $pair]]
}
puts $file [::csv::join $bondl]
unset bondl
for {set pair 16} {$pair < 32} {incr pair} {
     lappend bondl [measure bond [lindex $Pairs $pair]]
}
puts $file [::csv::join $bondl]
close $file
}
