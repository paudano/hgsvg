#!/usr/bin/env python

import sys
frac = float(sys.argv[1])
op =  sys.argv[2]
for line in sys.stdin:
    if line[0] == "#":
        sys.stdout.write(line)
        continue

    v = line.split()
    if op == "gte":
        if float(v[-1]) >= frac:
            sys.stdout.write(line)
    elif op == "lt":
        if float(v[-1]) < frac:
            sys.stdout.write(line)
