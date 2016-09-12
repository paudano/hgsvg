#!/usr/bin/env python

import argparse
import pysam
import Overlap
import Paths
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq

ap = argparse.ArgumentParser(description="Patch paths ")
ap.add_argument("ovp", help="overlaps")
ap.add_argument("asm", help="assemblies")
ap.add_argument("paths", help="Paths to patch")
ap.add_argument("contigs", help="Write contigs here.")
ap.add_argument("--chr", help="Process this chromosome.")
ap.add_argument("--patchbed", help="Patch contigs into chromosome, rather than separate ones using this bed file.", default=None)
ap.add_argument("--junction-N", help="Swap out this many bases at the junction between two contigs with N's. This can help debug joining, since misjoins will show large indels at the patch sites.",type=int, default=0, dest="junctionN")
                
args = ap.parse_args()


class SamBed:
    def __init__(self, line):
        v = line.split()
        self.chrom = v[0]
        self.tStart = int(v[1])
        self.tEnd   = int(v[2])
        self.read   = v[3]
        self.strand = int(v[4])
        self.qStart = int(v[5])
        self.qEnd   = int(v[6])

def ReadName(line):
    return line.split()[3]
        
def ReadSamBed(filename):
    samBedFile = open(filename)
    return { ReadName(line): SamBed(line) for line in samBedFile }


def ReplaceN(s,l):
    return s[-l:] +'N'*l


def PatchPath(ovpQuery, asm, path):
    #
    # Grab the first file.
    #
    if len(path) < 2:
        return
    edge = (path[0], path[1])
    segments = []
    if edge not in ovpQuery:
        print "ERROR, did not find an overlap in the path at 0"
        sys.exit(1)
        
    ovp = ovpQuery[edge]
    seqStart = ovp.aRead[0]
    seqEnd   = ovp.aOvp[1]
    seq = asm.fetch(reference=ovp.a, start=seqStart, end=seqEnd)
    # patch with N's to learn where the gap is 
    seq = ReplaceN(seq, args.junctionN)
    segments.append(seq)
    for i in range(1, len(path)-1):
        prevOverlapEdge = (path[i-1], path[i])
        if prevOverlapEdge not in ovpQuery:
            print "ERROR, did not find an overlap that is in the path at " + str(i)
            sys.exit(0)
            
        prevOverlap = ovpQuery[prevOverlapEdge]
        segmentStart = prevOverlap.bOvp[1]
        curOverlapEdge = (path[i], path[i+1])
        curOverlap = ovpQuery[curOverlapEdge]
        segmentEnd = curOverlap.aOvp[1]
        print "start: " + str(segmentStart) + " end: " + str(segmentEnd)
        if segmentEnd < segmentStart:
            import pdb
            pdb.set_trace()
        segment = asm.fetch(reference=curOverlap.a, start = segmentStart, end = segmentEnd)
        segment = ReplaceN(segment, args.junctionN)
        segments.append(segment)

    if len(path) > 2:
        edge = (path[-2], path[-1])
        
        ovp = ovpQuery[edge]
        segmentStart = ovp.bOvp[0]
        segmentEnd   = asm.get_reference_length(ovp.b) - ovp.bRead[1]
        segment = asm.fetch(reference=curOverlap.a, start = segmentStart, end = segmentEnd)
        segment = ReplaceN(segment, args.junctionN)
        segments.append(segment)

    contig = ''.join(segments)
    return seq


        

overlaps = Overlap.ReadOverlapFile(args.ovp)
asm      = pysam.FastaFile(args.asm)
paths    = Paths.ReadPaths(args.paths)

ovpQuery = Overlap.MakeOverlapQuery(overlaps)

contigOut = open(args.contigs,'w')

for pathIndex in range(0,len(paths)):
   contig = PatchPath(ovpQuery, asm, paths[pathIndex])
   if contig is None:
       continue
   contigSeq = SeqRecord.SeqRecord(Seq.Seq(contig), id="stitch.{:0>3}.{}.{}".format(pathIndex, len(paths[pathIndex]), len(contig)), name="", description="")

   SeqIO.write(contigSeq, contigOut, "fasta")

   
    
