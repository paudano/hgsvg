#!/usr/bin/env python
import sys
#1	chr1
#2	819573
#3	821888
#4	chr1
#5	820573
#6	820888
#7	DEL
#8	315
#9	NONE
#10	0,0
#11	0,0
#12	820888
#13	DEL
#14	11
#15	Pindel
#16	chr1
#17	820900
#18	820970
#19	deletion
#20	70
#21	chr1:796179-876179/0
#22	27301
#23	27301
#24	LOCAL
#25	HAP0
#


import argparse

ap = argparse.ArgumentParser(description="Select best sv overlap")
ap.add_argument("--infile", help="Input file",default="/dev/stdin")
ap.add_argument("--count", help="Count values with at least this overlap", default=None,type=float)
ap.add_argument("--header", help="Write header", default=None)
ap.add_argument("--bnidx", help="Use indexing for bionano", default=False, action="store_true")
ap.add_argument("--ts", help="Specify target start explicitly", default=None, type=int)
ap.add_argument("--te", help="Specify target end explicitly", default=None, type=int)
ap.add_argument("--qs", help="Specify query start",default=None,type=int)
ap.add_argument("--qe", help="Specify query end",default=None,type=int)
ap.add_argument("--svlen", help="Specify index of svlen", default=7,type=int)
ap.add_argument("--bntarget", help="Use bionano target", default=False, action="store_true")
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--key", help="Add key for best overlap", default=False, action="store_true")
ap.add_argument("--querykey", help="Add query key", default=None)
ap.add_argument("--operation", help="Filter raw reads by operation.", default=None)
ap.add_argument("--minOverlap", help="Minimum overlap", default=0.5, type=float)
args = ap.parse_args()

prev=None
outFile = open(args.out,'w')
if args.header:
    header="{}".format(args.header)
else:
    header="best"
if args.count:
    header+="\tnumber"
if args.key:
    header+= "\t{}".format(args.header +"_key")

if args.querykey:
    header += "\t" + args.querykey

if args.header is not None:
    outFile.write(header+"\n")

def Overlap(a, b):
    if b[0] == "-1" or b[0] == ".":
        return 0
    if args.bnidx or args.bntarget:
        # Compute bionano overlap -- just the overlap of lengths of events.
        
        al = float(int(a[1]) - int(a[0]))
        bl = float(int(b[1]) - int(b[0]))

        if min(al,bl) == 0:
            return 0
        else:
            return min(al,bl)/max(al,bl)
    
    if b[0] == ".":
        return 0
    a[0], a[1] = int(a[0]), int(a[1])
    b[0], b[1] = int(b[0]), int(b[1])
    if a[1] < b[0] or a[0] > b[1]:
        return 0
    isect = min(a[1],b[1]) - max(a[0],b[0])
    al=a[1]-a[0]
    bl=b[1]-b[0]
    if al > 0 and bl > 0:
        return min(float(al)/bl, float(bl)/al)
    
    return isect/float((max(b[1]-b[0], a[1]-a[0])))

matches=[]
infile = open(args.infile)
qs=4
qe=5
if args.qs is not None:
    qs = args.qs
if args.qe is not None:
    qe = args.qe
    
ts=16
te=17

    
if args.bnidx:
    qs=4
    qe=5
    ts=14
    te=15

if args.ts is not None:
    ts = args.ts
if args.te is not None:
    te = args.te


def GetEnd(m):
    if args.bnidx:
        if m[args.svlen] == "-1":
            return -1
        else:
            return int(m[qs]) + int(m[args.svlen])
    else:
        if m[qe] == "-1":
            return "-1"
        else:
            return m[qe]

def GetTargetEnd(m):
    if args.bntarget:
        if m[te+2] == "-1":
            return -1
        else:
            if int(m[te]) == (-1):
                return -1
            else:
                return int(m[ts]) + int(m[te+2])
    else:
        if int(m[te]) == -1:
            return -1
        else:
            return int(m[te])

def ProcessMatches(matches):
    if len(matches) == 0:
        res="0"
        if args.count:
            res+="\t0"
        if args.key:
            res+="\tNONE"
        if args.querykey:
            res+="\tNONE"
        outFile.write(res+"\n")
        return
    optMatch=[]
    filtMatches =[]
    for m in matches:
        if len(m) > te and m[ts] != "." :
            filtMatches.append(m)
    if len(filtMatches) == 0:
        res="0"
        if args.count:
            res+="\t0"
        
        if args.key:
            res+="\tNONE"
        if args.querykey:
            res+="\tNONE"
            
        outFile.write(res+"\n")
        
        return
    
    ovps = [ Overlap([m[qs], GetEnd(m)], [m[ts], GetTargetEnd(m)])  for m in filtMatches ]

    maxOverlap = max(ovps)
    maxOverlapIndex = 0
    nPass = 0
    if args.count is not None:
        for i in range(0,len(filtMatches)):
            if ovps[i] >= args.count:
                nPass+=1
    else:
        nPass = 1
    key="NONE"
    queryKey ="NONE"
    if maxOverlap > args.minOverlap:
        for i in range(0,len(filtMatches)):
            if ovps[i] == maxOverlap:
                maxOverlapIndex = i
                key = matches[i][ts-1] + "_" + matches[i][ts] + "_" + matches[i][te]
                queryKey = matches[i][0] + "_" + matches[i][qs] + "_" + matches[i][qe]
                break
            maxMatch = filtMatches[maxOverlapIndex]
    res = "{:2.2f}".format(maxOverlap)
    if args.count:
        res+="\t" + str(nPass)
    if args.key:
        res+= "\t"+key
    if args.querykey:
        res+="\t" + queryKey
        
    outFile.write(res+"\n")


def GetKey(vals, headerKey):
    key = "_".join(vals[0:3])
    if "hap" in headerKey:
        key+= "_"+vals[headerKey["hap"]]
    if "svLen" in headerKey:
        key+= "_"+vals[headerKey["svLen"]]
    return key

for line in infile:
    if line[0] == "#":
        headerVals = line.split()
        headerKey  = {headerVals[i]:i for i in range(0,len(headerVals))}
        prev=None
        if "oStart" in headerKey:
            ts = headerKey["oStart"]
        if "oEnd" in headerKey:
            te = headerKey["oEnd"]
            args.svlen = headerKey["svLen"]
        if "tStart" in headerKey:
            qs = headerKey["tStart"]
        if "tEnd" in headerKey:
            qe = headerKey["tEnd"]
        continue
    vals = line.split()
    key  = GetKey(vals, headerKey)

    if prev is not None and key != prev:
        ProcessMatches(matches)
        matches = []
    matches.append(vals)
    prev = key

#print headerKey


ProcessMatches(matches)
            
            
        
        
            
            
            
