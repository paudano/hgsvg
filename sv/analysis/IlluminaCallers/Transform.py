#!/usr/bin/env python
import sys
op = sys.argv[1]
for line in sys.stdin:
    if line[0] != "#" and op == "INS":
        vals = line.split()
        vals[2] = str(int(vals[1])+1)
        line = "\t".join(vals) + "\n"
    sys.stdout.write(line)
