#!/usr/bin/env python

import sys
def swap(a, i, j):
    t = a[i]
    a[i] = a[j]
    a[j] = t

for line in sys.stdin:
    vals = line.split()
    if len(vals) < 10:
        continue
    swap(vals, 0, 7)
    swap(vals, 1, 8)
    swap(vals, 2, 9)
    sys.stdout.write("\t".join(vals) + "\n") 
