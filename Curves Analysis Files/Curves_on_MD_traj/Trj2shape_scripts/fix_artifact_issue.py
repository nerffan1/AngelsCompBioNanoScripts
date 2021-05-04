import os 
import sys 
from collections import defaultdict 
from optparse import OptionParser 

parser = OptionParser() 
parser.add_option("-i", "--in", dest="input_filename",
                  help="input shape file ", metavar="FILE")
parser.add_option("-a", "--artifact_file", dest="artifact_filename",
                  help="artifact file, this file is generated why getting"
                       "shape feature values from lis files, check for"
                       " <.artifact> extension file in the directory where"
                       " you have shape files", metavar="FILE")
parser.add_option("-o", "--output_file", dest="output_filename",
                  help="output file. This is the file you should be"
                       "using for the downstream analysis", metavar="FILE")

(options, args) = parser.parse_args()
if len(sys.argv) == 1:
    parser.print_help()
    exit(1)

input_filename = options.input_filename 
artifact_filename = options.artifact_filename 
output_filename = options.output_filename  

artifact_files = [ l.strip() for l in open (artifact_filename) ]

input_file_lines = [l.strip() for l in open (input_filename) ]  

ofp = open (output_filename, "w") 

for l in input_file_lines: 
    lis_filename = l.split()[0] 
    if  lis_filename not in artifact_files: 
        ofp.write(l + "\n")
