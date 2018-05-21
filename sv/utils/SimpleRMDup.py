#!/usr/bin/env python
import sys
prev=None

def Ovp(a, b):
    ac= a[0]
    bc= b[0]
    ast =int(a[1])
    bs = int(b[1])
    ae = int(a[2])
    be = int(b[2])
    maxs=max(ast,bs)
    mins=min(ast,bs)
    maxe=max(ae,be)
    mine=min(ae,be)

    if ac != bc:
        return False

    if (maxs <= mine):
        sys.stderr.write("OVP!\n")
        return True
    return False
prev=None
for line in sys.stdin:
    if line[0] == "#":
        sys.stdout.write(line)
        continue

    if prev is None:
        prev = line.rstrip().split()
        continue
    vals = line.split()
    if Ovp(prev, vals):
        continue
    else:
        sys.stdout.write("\t".join(prev)+ "\n")
        prev=vals
if prev is not None:
    sys.stdout.write("\t".join(prev)+ "\n")        
