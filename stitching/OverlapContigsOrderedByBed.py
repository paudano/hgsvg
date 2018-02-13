#!/usr/bin/env python

import argparse
from multiprocessing import Process, Lock, Semaphore, Pool
import sys
import subprocess


ap = argparse.ArgumentParser(description="Given an alignmetn bed file, run overlaps between contigs that are near each other")
ap.add_argument("bed", help="Alignment bed file.")
ap.add_argument("asm", help="Assemblies file")
ap.add_argument("--chrom", help="Process this chromosome only.", default =None)
ap.add_argument("--path", help="Path to HGSVG stitching scripts.",\
                default="/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/stitching")
ap.add_argument("--ahead", help="Allow overlap betewen alignemnts ending this far ahead of current.",\
                type=int, default=20000)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--nproc", help="Number of processors.", default=1,type=int)
ap.add_argument("--tmpdir", help="Pass along tempdir.", default=".")
ap.add_argument("--blasr", help="Blasr command to run.", default=None)
ap.add_argument("--minOvp", help="Minimum overlap to store", default=10000,type=int)
ap.add_argument("--start", help="First line to start overlapping", default=0,type=int)
ap.add_argument("--n", help="N entries to compute overlap for", default=-1,type=int)
args = ap.parse_args()

def ParseBedLine(line):
    v = line.split()
    return [v[0], int(v[1]), int(v[2]), v[3]]

bedFile = open(args.bed)
outFile = open(args.out,'w')
alignments= [ParseBedLine(line) for line in bedFile]


if (args.n == -1):
    last = len(alignments)
else:
    last = min(len(alignments), args.start+args.n)
start=args.start
stepBack=start

while stepBack > 0 and \
      alignments[stepBack-1][0] == alignments[start][0]  and \
      alignments[start][1] <= alignments[stepBack][2]+args.ahead:
    stepBack-=1

start=stepBack

alignments = alignments[start:last]

commands = []
def Run(command):
    res = subprocess.check_output(command, shell=True)
    
    return res.rstrip()

if args.blasr is not None:
    blasrCommand = " --blasr " + args.blasr
else:
    blasrCommand = ""


n=0
firstTargetOffset=args.start-start
for i in range(1+firstTargetOffset,len(alignments)):
    if args.chrom is not None and alignments[i][0] != args.chrom:
        continue
    j = i-1
    while j > 0 < len(alignments) and \
          alignments[j][0] == alignments[i][0] and \
          alignments[j][2] + args.ahead > alignments[i][1]:
        j-=1

    aAln = " ".join([alignments[k][3] for k in range(j,i)])
    commands.append( "{}/OverlapTwoContigs.py --a \"{}\" --b {} --asm {} --tmpdir {} {} --minOvp {}".format(args.path, aAln, alignments[i][3], args.asm, args.tmpdir, blasrCommand, args.minOvp))

    
    n+=1
if len(commands) == 0:
    outFile.close()
    sys.exit(0)
#
# Now process alignments
#
pool = Pool(args.nproc)
ovps=[]
#res=pool.map_async(Run, commands, callback=ovps.append )
ovps = [ Run(c) for c in commands ]

filt=[]
def GetStarts(l):
    vals=l.split()
    v1=vals[0].split(":")[1].split("-")[0]
    v2=vals[3].split(":")[1].split("-")[0]
    return ((int(v1), int(v2)))


def Compare(a,b):
    sa=GetStarts(a)
    sb=GetStarts(b)
    if (sa[0] != sb[0]):
        return sa[0] < sb[0]
    else:
        return sa[1] < sb[1]

for ovp in ovps:
    if ovp != "\n" and ovp != "":
        filt.append(ovp)
    
filt.sort(cmp=Compare)

outFile.write("\n".join(filt)+"\n")
