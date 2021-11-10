# Description: This short script exchanges lines for a "standard_b.lib" file
# necessary for Curves+ input. It exchanges (using python list indexing), from
# 10 lines, 2 <--> 3, and then 2 <--> 9. But it needs to this for all sets of 10
# lines.
# Author: Angel-Emilio Villegas Sanchez
#Last Updated: 07/05/21

import sys


if len(sys.argv) == 2 :
    r_file = open(sys.argv[1],'r')
    lines = list(line.strip() for line in r_file)


    # Exchange the lines. Not necessary to alocate extra memory, so I'll just have
    # one counter.
    for sec in range(2,len(str)-1,10):
        third = sec+1
        lines[sec],lines[third] = lines[third],lines[sec]
        ninth = sec + 7
        lines[sec],lines[ninth] = lines[ninth],lines[sec]

    for lin in lines:
        print(lin)
