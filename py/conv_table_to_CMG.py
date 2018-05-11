"""conv_table_to_CMG - write CMG-friendly file of coordinates."""

import sys

outfile = str(sys.argv[1]) + '_Parsed.txt'
infile = str(sys.argv[1]) + '.txt'

g = open(outfile, 'w')

with open(infile, 'r') as f:
    next(f)
    for line in f:
        x, y, z = line.split()
        out = str(x)+' '+str(y)+' '+str(z)+' = 0\n'
        g.write(out)

g.close()
