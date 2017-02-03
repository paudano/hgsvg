#!/usr/bin/env python
import sys

sample = sys.argv[1]

for line in sys.stdin:
    if line[0:3] == "@RG":
        vals = line.split()
        vals[3] = "SM:" + sample
        sys.stdout.write("\t".join(vals) + "\n")
    else:
        sys.stdout.write(line)
