#!/usr/bin/env python

import sys
first = True
for line in sys.stdin:
    if first and line[0] == '#':
        vals = line.split()
        counts = {}
        for i in range(0,len(vals)):
            if vals[i] in counts:
                counts[vals[i]]+=1
                vals[i]+= "_"+str(counts[vals[i]])
            else:
                counts[vals[i]]=0

        sys.stdout.write("\t".join(vals) + "\n")
        first = False
    else:
        sys.stdout.write(line)
                
