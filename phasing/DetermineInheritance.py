#!/usr/bin/env python

import pysam
import argparse

ap = argparse.ArgumentParser(description="Find which haplotype was inherited, for a region")
ap.add_argument("--vcf", help="VCF file, should have trio", required=True)
ap.add_argument("--fa", help="Sample ID of father to consider.", required=True)
ap.add_argument("--mo", help="Sample ID of mother to consider.", required=True)
ap.add_argument("--child", help="Sample ID of child", required=True)
#ap.add_argument("--region", help="Region to seek phase.", required=True)
ap.add_argument("--faBed", help="Output file.", required=True)
ap.add_argument("--moBed", help="Output file.", required=True)
args = ap.parse_args()


vcfFile = pysam.VariantFile(args.vcf)  # auto-detect input format


#rgn = args.region.split(".")
faOut = open(args.faBed, 'w')
moOut = open(args.moBed, 'w')
def GetInheritance(gt, ph):
    if ph[0] == gt:
        return (0,1)
    return (1,0)

faStart = None
moStart = None
moPrev  = None
faPrev = None
for rec in vcfFile.fetch():
    #
    # Father is homozygous, motheris het, so it is possible to determine hertiance from mother
    #

    if rec.samples[args.child]['PS'][0] != '.' and \
        (rec.samples[args.mo]['PS'][0] != '.' and \
         rec.samples[args.fa]['PS'][0] == '.'):
         ch = rec.samples[args.child]['GT']
         fa = rec.samples[args.fa]['GT']
         mo = rec.samples[args.mo]['GT']
         if fa[0] == fa[1]:
             faGT = fa[0]
             (faHap,moHap) = GetInheritance(faGT, ch)
             (moInherited, moRetained) = GetInheritance(ch[moHap], mo)
             if moPrev is None or moPrev[0] != rec.chrom:
                 if moPrev is not None:
                     moOut.write(str(moStart[0]) + "\t" + str(moStart[1]) + "\t" + str(moPrev[1]) +"\t" + str(moPrev[2]) + "\t" + str(moPrev[3]) + "\n")                     
                 moPrev = (rec.chrom, rec.start, moInherited, moHap)
                 moStart = (rec.chrom, rec.start, moInherited, moHap)
             else:
                 moCur = (rec.chrom, rec.start, moInherited, moHap)
                 if moCur[2] != moPrev[2]:
                     moOut.write(str(rec.chrom) + "\t" + str(moStart[1]) + "\t" + str(moPrev[1]) +"\t" + str(moPrev[2]) + "\t" + str(moPrev[3]) + "\n")
                     moStart = moCur
                 moPrev = moCur
                 
    elif rec.samples[args.child]['PS'][0] != '.' and \
         rec.samples[args.mo]['PS'][0] == '.' and \
         rec.samples[args.fa]['PS'][0] != '.':
         ch = rec.samples[args.child]['GT']
         fa = rec.samples[args.fa]['GT']
         mo = rec.samples[args.mo]['GT']
         if mo[0] == mo[1]:
             moGT = mo[0]
             (matChrom,patChrom) = GetInheritance(moGT, ch)
             (faPassed, faRetained) = GetInheritance(ch[patChrom], fa)

             if faPrev is None or faPrev[0] != rec.chrom:
                 if faPrev is not None:
                     faOut.write(str(faStart[0]) + "\t" + str(faStart[1]) + "\t" + str(faPrev[1]) +"\t" + str(faPrev[2]) + "\t" + str(faPrev[3]) + "\n")
                 faPrev = (rec.chrom, rec.start, faPassed, patChrom)
                 faStart = (rec.chrom, rec.start, faPassed, patChrom)
             else:
                 faCur = (rec.chrom, rec.start, faPassed, patChrom)
                 if faCur[2] != faPrev[2]:
                     faOut.write(str(rec.chrom) + "\t" + \
                                 str(faStart[1]) + "\t" + \
                                 str(faPrev[1]) +"\t" + \
                                 str(faPrev[2]) + "\t" + \
                                 str(faPrev[3]) + "\n")
                     faStart = faCur
                 faPrev = faCur

moOut.write(str(rec.chrom) + "\t" + str(moStart[1]) + "\t" + str(moPrev[1]) + "\t" + str(moPrev[2])+ "\t" + str(moPrev[3]) + "\n")
faOut.write(str(rec.chrom) + "\t" + str(faStart[1]) + "\t" + str(faPrev[1]) + "\t" + str(faPrev[2])+ "\t" + str(faPrev[3]) + "\n")
        
faOut.close()
moOut.close()
