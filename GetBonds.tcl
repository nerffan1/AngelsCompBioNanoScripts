# GET GetGGBonds
#------------------------------------------------------------------
# Description: The GetGGBonds procedure gets the bonds from all the pdbs in a directory.
# Each bond will be appended to a csv file on a line

package require csv

proc GetGGBonds {} {
     set first 1
     foreach struct [glob *.pdb] {
          mol load pdb $struct
          set sel [atomselect top "all"]
          if {$first == 1} {
               set Pairs [GetBondPairs $sel]
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

#Description: After obtaining a list with all pair of indices corresponding to bonds, WriteCSV
# appends to a file the bond values
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
