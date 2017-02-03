#!/usr/bin/env python

import argparse
from Bio import Seq
from Bio import SeqRecord
from Bio import SeqIO
import pysam
import tempfile
import subprocess
from multiprocessing import Process, Lock, Semaphore, Pool

import re
import os
import sys

ap = argparse.ArgumentParser(description="Given a gaps.bed file, splice variants into the reference, remap reads to this, and count gaps. ")
ap.add_argument("--gaps", help="Gaps file (can be insertions and deletions).", required=True)
ap.add_argument("--ref", help="Reference file.", required=True)
ap.add_argument("--reads", help="FOFN of aligned reads.", required=True)
ap.add_argument("--flank", help="Use this amount of the reference when aligning reads", default=900,type=int)
ap.add_argument("--window", help="Collect all gaps within this window.", type=int, default=1000)
ap.add_argument("--minGap", help="Minimum gap size.", default=30, type=int)
ap.add_argument("--out", help="Output file.", default="/dev/stdout")
ap.add_argument("--tmpdir",help="Temporary directory",default=None)
ap.add_argument("--nproc", help="Number of threads.", default=4,type=int)
ap.add_argument("--keep", help="Keep temporary files", default=False, action='store_true')
ap.add_argument("--mssm", help="Parse mssm lines.", default=False, action='store_true')
ap.add_argument("--falcon", help="Parse in falcon format", default=False, action='store_true')
ap.add_argument("--count", help="Just count clusters, don't do any validation.", default=None)
ap.add_argument("--haps", help="Consider these haplotypes", default=[], nargs="+")

args = ap.parse_args()
if args.count is None:
    outFile =open(args.out,'w')

if args.tmpdir is None:

    if "TMPDIR" not in os.environ or os.environ["TMPDIR"] == "":
        print "ERROR. The TEMPDIR variable must be set or --tmpdir specified on as a command  argument"
        sys.exit(1)
    else:
        args.tmpdir = os.environ["TMPDIR"]

ref = pysam.FastaFile(args.ref)

refFaiFile = open(args.ref + ".fai")
refFai = {l.split()[0]: int(l.split()[1]) for l in refFaiFile}
readsFofn = open(args.reads)
bamFiles = [pysam.AlignmentFile(line.rstrip(), 'rb') for line in readsFofn]

def CountRefCoverage(samFileName):
    samFile = pysam.AlignmentFile(samFileName,'r')
    refs = {}
    refs={}
    refs["re"] = 0
    refs["db"] = 0
    for ref in samFile.references:
        refs[ref] = 0

    for read in samFile.fetch():
        if read.reference_id != -1:
            rn = read.reference_name 
            if rn not in refs:
                refs[rn] = 0
            refs[rn]+=1
    return refs

def CountGaps(samFileName, rStart, rEnd):
    samFile = pysam.AlignmentFile(samFileName,'r')
    nGaps = 0
    nReads = 0
    nGapBases=0
    for read in samFile.fetch():
        nReads+=1
        
        if read.cigartuples is not None:
            refPos = 0
            for t in read.cigartuples:
                
                if (rStart  <= refPos and refPos < rEnd):
                    if (t[0] == 1 or t[0] == 2) and t[1] > args.minGap:
                        nGaps+=1
                        nGapBases+=t[1]
                    if t[0] == 4:
                        nGaps+=1

                if t[0] == 0 or t[0] == 2 or t[0] == 7 or t[0] == 8:
                    refPos+=t[1]
    return nGaps, nGapBases, nReads


def CountCoverage(samFileName, rStart, rEnd):
    samFile = pysam.AlignmentFile(samFileName,'r')
    nGaps = 0
    nReads = 0
    nCovered = 0
    for read in samFile.fetch():
        nReads+=1
        covered=False
        hasGap =False
        
        if read.cigartuples is not None:
            refPos = 0
            for t in read.cigartuples:
                
                if (rStart  <= refPos and refPos < rEnd):
                    if (t[0] == 1 or t[0] == 2) and t[1] > args.minGap:
                        hasGap=True
                        break
                    if t[0] == 4:
                        hasGap=True
                        break

                if t[0] == 0 or t[0] == 2 or t[0] == 7 or t[0] == 8:
                    refPos+=t[1]
        if hasGap == False:
            nCovered+=1
    return nCovered, nReads


class Gap:
    def __init__(self):
        self.chrom=""
        self.start= 0
        self.end  = 0
        self.svType=0
        self.svLen=0
        self.hap=None
        self.rem=[]

def ParseFalconLine(vals):
    gap=Gap()
    gap.chrom=vals[0]
    gap.start=int(vals[1])
    gap.end=int(vals[2])
    gap.svType="deletion"
    return gap

def ParseMSSMLine(vals):
    print "parsing mssm!"
    if len(vals) > 0:
        gap=Gap()
        gap.chrom=vals[0]
        gap.start=int(vals[1])
        gap.end=int(vals[2])
        gap.svType=vals[3]
        gap.svLen=len(vals[5])
        return gap
    else:
        return None

def ParseGapsLine(vals):
    if len(header) > 0:
        chromI = header["chrom"]
        tStartI = header["tStart"]
        tEndI   = header["tEnd"]
        svTypeI = header["svType"]
        svLenI   = header["svLen"]
        svSeqI = header["svSeq"]
        remainder = 6
    else:
        chromI, tStartI, tEndI, svTypeI, svLenI, svSeq = 0,1,2,3,4,5

        
    if len(vals) > 0:
        gap=Gap()
        gap.chrom=vals[chromI]
        gap.start=int(vals[tStartI])
        gap.end=int(vals[tEndI])
        gap.svType=vals[svTypeI]
        gap.seq=vals[svSeqI]

        gap.svLen=int(vals[svLenI])
        if "hap" in header:
            gap.hap=vals[header["hap"]]
        else:
            gap.hap=None
        return gap

    else:
        return None


def GetEnd(v):
    if args.mssm:
        print v.svLen
        return v.start + v.svLen
    elif args.falcon:
        return v.end
    else:
        if v.svType == "deletion":
            return v.end
        else:
            return v.start+1

def GetStart(v):
    return v.start

def GetChrom(v):
    return v.chrom

gapsFile = open(args.gaps)
headerLine = gapsFile.readline()
gapLines = gapsFile.readlines()

if headerLine[0] != '#':
    gapLines = [header] + gapLines
    header=None
else:
    v = headerLine[1:].split()
    header = {v[i]: i for i in range(0,len(v))}

if args.mssm:
    gaps = [ParseMSSMLine(line.split()) for line in gapLines]
elif args.falcon:
    gaps = [ParseFalconLine(line.split()) for line in gapLines]
else:
    gaps = [ParseGapsLine(line.split()) for line in gapLines]

temp=[]
for g in gaps:
    if g is not None:
        temp.append(g)
gaps= temp

svClusters = []
i=0
while i < len(gaps):
    svClusters.append([])
    svClusters[-1].append(gaps[i])

    while i < len(gaps)-1 and GetChrom(gaps[i]) == GetChrom(gaps[i+1]) and GetStart(gaps[i+1])-GetEnd(gaps[i]) < args.window:
        i+=1
        svClusters[-1].append(gaps[i])
    i+=1

if args.count is not None:
    countFile = open(args.count,'w')
    countFile.write("\n".join([str(len(c)) for c in svClusters]))
    countFile.close()
    sys.exit(0)
    
def WriteSeq(fh,seq,seqName):
    SeqIO.write(SeqRecord.SeqRecord(Seq.Seq(seq),id=seqName,name="",description=""),fh,"fasta")

sem = Semaphore(1)

def SpliceTestLine(svs):
    svIndex = 0
    for sv in svs:
        if sv.chrom not in refFai:
            continue
        if svIndex == 0:
            prefixStart = max(0, sv.start-args.flank)
            prefixEnd   = sv.start
            sem.acquire()
            gapPrefix   = ref.fetch(sv.chrom, prefixStart, prefixEnd)
            sem.release()
            spliceSeqs  = [gapPrefix]
#            print "gap prefix " + gapPrefix
            refPos      = sv.start
        
        if args.mssm:
            spliceSeqs.append(sv[5])

            refPos = sv[2]

            if svIndex < len(svs)-1:
                if refPos > svs[svIndex+1][1]:
                    continue
                sem.acquire()
                between = ref.fetch(sv[0], refPos, svs[svIndex+1][1])
                sem.release()
                spliceSeqs.append(between)
        elif args.falcon:
            # only checking deletions
            refPos = sv[2]
            if svIndex < len(svs)-1:
                if refPos > svs[svIndex+1][1]:
                    continue
            
        else:
            if sv.svType == "insertion":
#                print "appending: " + sv[5]
                spliceSeqs.append(sv.seq)
                refPos = sv.start
            else:
                refPos = sv.end
        

            if svIndex < len(svs)-1:
                if refPos > svs[svIndex+1].start:
                    sys.stderr.write("Ignoring an sv.\n")
                    sys.stderr.write(str(svs) + "\n")
                    continue
                else:
                    between = ref.fetch(sv.chrom, refPos, svs[svIndex+1].start)
#                    print "between " + between
                    refPos += len(between)
                    spliceSeqs.append(between)
        svIndex +=1

    svIndex = len(svs)-1
    sem.acquire()
    suffix = ref.fetch(svs[svIndex].chrom, refPos, min(refFai[svs[svIndex].chrom], refPos + args.flank))
    sem.release()
    
#    print "gap suffix " + suffix
    spliceSeqs.append(suffix)
    dbSeq= ''.join(spliceSeqs)
#    print str(len(dbSeq))
    refChrom = svs[0].chrom


    refStart = max(0, svs[0].start - args.flank)
    refEnd   = min(refFai[svs[svIndex].chrom], GetEnd(svs[svIndex]) + args.flank)
#    print "prefix: " + str(prefixStart) + " " + str(prefixEnd) + " " + str(refStart) + " " + str(refEnd)
    sem.acquire()
    refSeq = ref.fetch(refChrom, refStart, refEnd)
    sem.release()
    rFile     = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".fasta", delete=False, mode='w')
    dbFile    = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".fasta", delete=False, mode='w')
    readsFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".fasta", delete=False, mode='w')

    
    WriteSeq(dbFile, dbSeq, "db")
    WriteSeq(dbFile, refSeq, "re")
    WriteSeq(rFile, refSeq, "re")    
    #
    # Now collect all of the sequences.
    #
    fetchStart = svs[0].start - args.flank
    fetchEnd   = GetEnd(svs[-1]) + args.flank
#    print "Fetching " + svs[0][0] + ":" + str(fetchStart) + "-" + str(fetchEnd) + " ref: " + refChrom + ":" + str(refStart) + "-" + str(refEnd) +  " " + str(svs[0][1])

    sem.acquire()
    for b in range(0,len(bamFiles)):
        for read in bamFiles[b].fetch(sv.chrom, fetchStart, fetchEnd+1):
            WriteSeq(readsFile, read.seq, read.query_name)
    sem.release()
    readsFile.close()
    dbFile.close()
    rFile.close()
#    rsFile  = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".sam", delete=False, mode='w')
    dbsFile = tempfile.NamedTemporaryFile(dir=args.tmpdir, suffix=".sam", delete=False, mode='w')
    
    commandOptions = " -sam -bestn 1 -affineOpen 5 -affineExtend 5 -nproc 6 -minAlignLength {} ".format(int(1.5*args.flank))

    dbCommand = "/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr {} {} -clipping soft -out {} ".format(readsFile.name, dbFile.name, dbsFile.name) + commandOptions

    dn = open(os.devnull)
#    subprocess.call(rCommand.split(),stderr=dn)
    subprocess.call(dbCommand.split(),stderr=dn)
    cov = CountRefCoverage(dbsFile.name)
    if "db" in cov and cov["db"] < 500:
        print "spliced " + str(len(svs))        
        print "cov:"
        print cov
        print "{}:{}-{}".format(svs[0].chrom,svs[0].start,svs[0].end)

#    rGaps, rGapBases, rReads = CountGaps(rsFile.name, args.flank - 100, len(refSeq) - args.flank+100)
#    print "db:"
#    dbGaps, dbGapBases, dbReads = CountGaps(dbsFile.name, args.flank - 100, len(dbSeq) - args.flank+100)

#    rCov,  rReads = CountCoverage(rsFile.name, args.flank - 100, len(refSeq) - args.flank+100)
#    print "db:"
#    dbCov, dbReads = CountCoverage(dbsFile.name, args.flank - 100, len(dbSeq) - args.flank+100)

    
#    sys.stderr.write("processed {}:{}-{}\n".format(sv[0],sv[1],sv[2]))
#    cleanup = "/bin/rm {} {} {} {} {}".format(readsFile.name, rFile.name, dbFile.name, rsFile.name, dbsFile.name)
    cleanup = "/bin/rm {} {} {} {}".format(readsFile.name, rFile.name, dbFile.name, dbsFile.name)


#    do = open(os.devout)
#    subprocess.call(chkref.split(), stdout=sys.stdout)
    if args.keep is False:
        subprocess.call(cleanup.split())
    else:
        print cleanup

    dbCov = cov["db"]
    rCov  = cov["re"]
    
    results="\n".join(["{}:{}-{}\t{}\t{}".format(sv.chrom,sv.start,sv.end,dbCov,rCov) for sv in svs])
    return results
if header is not None:
    outFile.write("#region\tnAlt\tnRef\n")
if args.nproc > 1:    
    pool = Pool(args.nproc)

    res = pool.map(SpliceTestLine, svClusters)
    pool.close()
    outFile.write("\n".join(res)+"\n")
else:
    for svCluster in svClusters:
        res = SpliceTestLine(svCluster)
        outFile.write(res+ "\n")
