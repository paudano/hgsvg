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
ap.add_argument("--ref", help="Reference genome", default="/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta")
ap.add_argument("--repack", help="Repack smaller indels", default=False, action="store_true")
ap.add_argument("--window", default=1000,type=int)
ap.add_argument("--op", help="operation", default=None)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--nproc", help="Number of threads.", default=1,type=int)
ap.add_argument("--header", help="Write header", default=False,action="store_true")
ap.add_argument("--isbn", help="Is Bionano query", default=False,action="store_true")
ap.add_argument("--bnop", help="Bionano op", default=None)
ap.add_argument("--keeptemp", help="Keep temporary files.", default=False, action="store_true")
ap.add_argument("--ts", help="Configure target start index", type=int,default=16)
ap.add_argument("--te", help="Configure target end index", type=int, default=17)
args = ap.parse_args()

tmpdir="."
if "TMPDIR" in os.environ:
	 tmpdir=os.environ["TMPDIR"]
callFile = open(args.calls)
readsFofn = open(args.fofn)
blines = readsFofn.readlines()
bamFiles = [pysam.AlignmentFile(bline.rstrip(), 'rb') for bline in blines]
bamFileIndex = {}
sem = Semaphore(1)
for b in range(0,len(blines)):
    vals = blines[b].split(".")
    for v in vals:
        if "chr" in v:
            bamFileIndex[v]= b
            break

counter = Value("d", 0)

outFile = open(args.out,'w')


if args.header:
    outFile.write("overlap\tnumber\n")

class Status:
    def __init__(self):
        self.lineIndex=0
        self.nLines = 0
status=Status()

def ProcessLine(line):
    if line[0] == "#":
        return "maxOverlap\tnReads"
    vals= line.split()
    region = vals[0] + ":" + vals[1] + "-" + vals[2]
    
    tempFileNames = []
    
    samFile = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".sam", delete=False, mode='w')
    tempFileNames.append(samFile.name)
    samFile.close()
    
    samHandle = pysam.AlignmentFile(samFile.name, 'wh', header=bamFiles[0].header)
    if int(vals[2]) - int(vals[1]) > 10000:
        sys.stderr.write("SKIPPING " + "\t".join(vals)  + "\n")
        return [0,0]


    if vals[0] not in bamFileIndex:
        sys.stderr.write("Missing " + vals[0] + "\n")
        return [0,0]
    b = bamFileIndex[vals[0]]

    sem.acquire()
    sys.stderr.write("Collecting alignments " + str(counter.value) + " of " + str(status.nLines) + "\n")
    counter.value+=1
    try:
        covall = bamFiles[b].count_coverage(vals[0], int(vals[1]), int(vals[2]), quality_threshold=0)
    except:
        sys.stderr.write("ERROR reading from bam file " + readsFofn[b] + "\n")
        sem.release()
        return [0,0]
            
    maxCov = max(covall[0])
    sys.stderr.write("MAX COVERAGE: " + str(maxCov) + "\n")
    if maxCov > 200:
        sem.release()
        return [0,0]
    
    for read in bamFiles[b].fetch(vals[0], int(vals[1]), int(vals[2])):
        samHandle.write(read)
    sem.release()
    			 
    							 
    querySam = samFile.name
    if args.repack:
        sys.stderr.write("repacking\n")
        smallPackedSam = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".sam", delete=False, mode='w')
        cmd = "/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/sv/utils/RawReadGaps.py --infile {}".format(samFile.name)
        tempFileNames.append(smallPackedSam.name)
        subprocess.call(cmd.split(), stdout=smallPackedSam)
        smallPackedSam.close()
        largePackedSam = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".sam", delete=False, mode='w')
        tempFileNames.append(largePackedSam.name)
        cmd = "/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/sv/utils/MergeGaps.py --infile {} --gap -50 --mmr 0.05".format(smallPackedSam.name)
        subprocess.call(cmd.split(), stdout=largePackedSam)
        querySam=largePackedSam.name
        largePackedSam.close()
    sys.stderr.write("Collecting gaps\n")
    gaps = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".bed", delete=False, mode='w')
    tempFileNames.append(gaps.name)
    cmd = "/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintGaps.py {} {} --condense 50 --outFile {} ".format(args.ref, querySam, gaps.name)
    sys.stderr.write(cmd+"\n") 
    subprocess.call(cmd.split())

    gapsFileName = gaps.name
    if args.op is not None:
        op = args.op
        if args.op == "automatic":
            op = vals[header["svType"]]

        opGaps = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".bed", delete=False, mode='w')
        cmd = 'egrep \"^#|{}\" {} > {}'.format(op, gaps.name, opGaps.name)
        sys.stderr.write(cmd + "\n")
        p=subprocess.Popen(cmd, shell=True)
        p.wait()
        opGaps.flush()
        opGaps.close()
        tempFileNames.append(opGaps.name)
        gapsFileName = opGaps.name
    
    sys.stderr.write("Intersection\n")
    query = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".bed", delete=False, mode='w')
    query.write("#queryChrom\tqueryStart\tqueryEnd\t" + headerLine+"\n")
    query.write("{}\t{}\t{}\t{}\n".format(vals[0], int(vals[1])-args.window, int(vals[2]) + args.window, line))
    query.close()
    tempFileNames.append(query.name)
    isect = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".bed", delete=False, mode='w')
    cmd = 'h=`head -1 {}`; echo -e $h"\t"oChrom"\t"oStart"\t"oEnd > {}'.format(query.name, isect.name)
    p=subprocess.Popen(cmd, shell=True)
    p.wait()
    
    cmd = "bedtools intersect -a {} -b {} -loj >> {}".format(query.name, gapsFileName, isect.name)
    sys.stderr.write(cmd + "\n")
    p=subprocess.Popen(cmd, shell=True)
    p.wait()
    #subprocess.call(cmd.split(), stdout=isect)
    tempFileNames.append(isect.name)
    isect.close()
    ovp = tempfile.NamedTemporaryFile(dir=tmpdir, suffix=".ovp", delete=False, mode='w')

    if args.bnop is not None:
        bnop = "--op " + args.bnop

    bnIdxOp = ""
    if args.isbn is True:
        bnIdxOp = " --bnidx"
        
    cmd = "/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/sv/analysis/IlluminaCallers/SelectBestcall.py --infile {} --count 0.5 --ts {} --te {} --svlen 10 {} --out {}".format(isect.name, args.ts, args.te, bnIdxOp, ovp.name)
    tempFileNames.append(ovp.name)
    sys.stderr.write(cmd + "\n")

    subprocess.call(cmd.split())
    ovp.close()
    ovpFile = open(ovp.name)
    lines = ovpFile.readlines()
    if len(lines) > 0:
        return lines[0].split()
    else:
        return [0,0]
    tempFileNames.append(ovp.name)
    if args.keeptemp is False:
        cleanup = "/bin/rm  -f " + " ".join(tempFileNames)
        subprocess.call(cleanup.split())

lines = callFile.readlines()
headerLine = ""
if len(lines) > 0 and lines[0][0] == "#":
    headerLine = lines[0]
    headerVals = lines[0][1:].rstrip().split()
    header = {headerVals[i] : i for i in range(0,len(headerVals)) }
    lines = lines[1:]

status.nLines = len(lines)
if args.nproc > 1:    
    pool = Pool(args.nproc)

    res = pool.map(ProcessLine, lines)
    pool.close()
    outFile.write("\n".join(["{}\t{}".format(r[0],r[1]) for r in res])+"\n")
else:
    for line in lines:
        res = ProcessLine(line)
        outFile.write(str(str(res[0]) + "\t" + str(res[1]) + "\n"))
