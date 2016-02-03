import sys

outfile = str(sys.argv[1]) + '_Parsed.txt'
infile = str(sys.argv[1]) + '.txt'

g = open(outfile,'w')

with open(infile,'r') as f:
 next(f)
 for line in f:
  x,y,z = line.split()
  out = str(x)+':'+str(x)+' '+str(y)+':'+str(y)+' '+str(z)+':'+str(z)+' = 0\n'
  g.write(out)
  
g.close()
