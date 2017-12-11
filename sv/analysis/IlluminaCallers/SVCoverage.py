#!/usr/bin/env python

import argparse
import tempfile
import os
import sys
import pysam
import subprocess
from multiprocessing import Process, Lock, Semaphore, Pool, Value

ap = argparse.ArgumentParser(description="Determine how many pb reads overlap an Illumina call")
ap.add_argument("--calls", help="Calls", required=True)
ap.add_argument("--fofn", help="FOFN of pb-bams", required=True)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--nproc", help="Number of threads.", default=1,type=int)
ap.add_argument("--op", help="Guide coverage by operation", default="DEL")
ap.add_argument("--bionano", help="Compute bionano estimage of coverage -- min coverage in window specified by svLen for deletions, and max for insertions", action="store_true", default="False")
ap.add_argument("--writecov", help="Write coverae to this file.",default=None)
ap.add_argument("--svlen", help="Specify column of svlen.",type=int,default=4)
ap.add_argument("--header", help="Write header", default=False,action="store_true")
args = ap.parse_args()

tmpdir="."
if "TMPDIR" in os.environ:
	 tmpdir=os.environ["TMPDIR"]
callFile = open(args.calls)
readsFofn = open(args.fofn)
blines = readsFofn.readlines()
bamFiles = [pysam.AlignmentFile(bline.rstrip(), 'rb') for bline in blines]
bamFileIndex = {}
for b in range(0,len(blines)):
    vals = blines[b].split(".")
    for v in vals:
        if "chr" in v:
            bamFileIndex[v]= b
            break

outFile = open(args.out,'w')
sem = Semaphore(1)

if args.header:
    outFile.write("coverage\n")

if args.writecov is not None:
    covFile = open(args.writecov, 'w')
counter = Value("d", 0)

def ProcessLine(line):
    if line[0] == "#":
        return 0
    
    vals= line.split()
    regionStart = int(vals[1])
    regionEnd   = int(vals[2])
    region = vals[0] + ":" + vals[1] + "-" + vals[2]
    
    tempFileNames = []
    regionLength = regionEnd-regionStart
    cov = [0]*(regionLength)

    sys.stderr.write("Getting coverage for " + str(counter.value) + " of " + str(nsv) + " " + str(regionEnd - regionStart) + "\n")
    counter.value += 1
    if vals[0] not in bamFileIndex:
        sys.stderr.write("Missing " + vals[0] + "\n")
        return 0
    b = bamFileIndex[vals[0]]

    sem.acquire()
    covall = bamFiles[b].count_coverage(vals[0], regionStart, regionEnd, quality_threshold=0)
    sem.release()
    for i in range(0,len(covall[0])):
        cov[i] = covall[0][i] + covall[1][i] + covall[2][i]  + covall[3][i]

    
    
    total=sum(cov)
    if args.writecov is not None:
        covFile.write("\t".join([str(c) for c in cov]) + "\n")
        covFile.write("\t".join([str(p) for p in range(regionStart,regionEnd)]) + "\n")
        
    if args.bionano:
        svLen = int(vals[args.svlen])
        
        if len(cov) < svLen:
            if len(cov) == 0:
                return 0
            else:
                return float(sum(cov))/len(cov)
        else:
            minWSum = 0
            maxWSum = 0

            curSum = sum(cov[0:svLen])
            minWSum = curSum
            maxWSum = curSum
            for i in range(1, len(cov) - svLen ):
#                wSum=sum(cov[i:i+svLen])

                curSum-=cov[i-1]
                curSum+=cov[i+svLen-1]
                if curSum < minWSum:
                    minWSum = curSum
                if curSum > maxWSum:
                    maxWSum = curSum
#            sys.stderr.write(" ".join(["{:2.2f}".format(i) for i in mc]) + "\n")
            if svLen > 0:
                sys.stderr.write(str(svLen) + " minsv: " + str(minWSum) + " total " + str(sum(cov)) + "\tmin {:2.2f}\tmax {:2.2f}\tregion {:2.2f}".format(minWSum / float(svLen), maxWSum/float(svLen), sum(cov) / float(len(cov))) + " " + str(len(cov))+ "\n")
            if args.op == "DEL":
                if svLen > 0:
                    return minWSum/float(svLen)
                else:
                    return 0
            elif args.op == "INS":
                if svLen > 0:
                    return maxWSum/float(svLen)
                else:
                    return 0
            else:
                print "ERROR WITH READ DEPTH"
                sys.exit(1)
            
    avg = 0
    if regionLength > 0:
        avg=float(total)/regionLength

    return avg

lines = callFile.readlines()

if len(lines) > 0 and lines[0][0] == "#":
    lines = lines[1:]
nsv = len(lines)
if args.nproc > 1:    
    pool = Pool(args.nproc)
    res = pool.map(ProcessLine, lines)
    pool.close()
    outFile.write("\n".join(["{:2.2f}".format(r) for r in res])+"\n")
else:
    for line in lines:
        res = ProcessLine(line)
        outFile.write("{:2.2f}\n".format(res))
