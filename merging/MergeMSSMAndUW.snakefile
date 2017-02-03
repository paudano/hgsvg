import os
import tempfile

#
# A little complicated to find the temp dir
#
SSD_TMP_DIR = "/data/scratch/ssd"
if "TMPDIR" in os.environ:
    TMPDIR = os.environ['TMPDIR']
elif "TMPDIR" in config:
    TMPDIR = config['TMPDIR']
elif os.path.exists(SSD_TMP_DIR):
    TMPDIR = SSD_TMP_DIR
else:
    TMPDIR = tempfile.gettempdir()

configfile: "merge.json"


SLOP_FOR_SV_SEQUENCE_POSITIONS = 5000


faiFile = open(config['ref']+".fai")
chroms = [l.split()[0].rstrip() for l in faiFile]

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
cwd=os.getcwd()

shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

haps=["h0","h1"]
parents=["fa", "mo"]

allMSSM = [config["mssm-h0"],config["mssm-h1"],config["mssm-dn"]]

allMSSMBase = {f: os.path.basename(f) for f in allMSSM}
mssmBaseToPath = {os.path.basename(f): f for f in allMSSM}
sources=["uw","mssm"]
longOps=["deletion","insertion"]
shortOps=["del","ins"]
filtered=["uw.bed.filt", "mssm.bed.filt"]

print(config["bionano"])
print(sources)
rule all:
    input:
        mssmSort   = expand("{mssm}.sorted",mssm=allMSSM),    
        mssmBPSupport=expand("{mssm}.cov",mssm=allMSSM),
        mssmToBN   = expand("{mssm}.to-bn.bed",mssm=mssmBaseToPath.keys()),
        mssmToUW   = expand("{mssm}.to-uw.bed",mssm=mssmBaseToPath.keys()),
        mssmToUWFilt = expand("{mssm}.to-uw.bed.filt",mssm=mssmBaseToPath.keys()),
        mssmBNAnnotation   = expand("{mssm}.to-uw.bed.bn",mssm=mssmBaseToPath.keys()),
        uwBNAnnotation = "uw.bed.bn",
        mergedMSSM = "mssm.bed.filt",
        uwFilt="uw.bed.filt",
        svbb=expand("{source}.{op}.bb", source=sources, op=longOps),
        svs=expand("{filt}.{op}.bed", filt=filtered,op=shortOps),
        svsbn=expand("{filt}.{op}.bed.bn", filt=filtered,op=shortOps),
        overlaps=expand("{filt}.{op}.overlaps.bed", filt=filtered,op=shortOps),
        svAnnot=["uw.bed.filt.ins.mssm.bed.filt-bn.bed","uw.bed.filt.del.mssm.bed.filt-bn.bed", "mssm.bed.filt.ins.uw.bed.filt-bn.bed","mssm.bed.filt.del.uw.bed.filt-bn.bed"],
        svDist=expand("mssm.bed.filt.{op}.bed.uw-dist",op=shortOps),
        svUWDist=expand("uw.bed.filt.{op}.bed.mssm-dist",op=shortOps),
        mergedCalls=expand("merged.{op}.bed",op=shortOps),
        mergedPreInvBed="sv_calls.pre-inv.bed",
        clusterCount=[config["uwsv"]+".cluster-count"] + expand("{mssm}.to-uw.bed.cluster-count",mssm=mssmBaseToPath.keys()),
        clusterSummary="cluster-summary.txt",
        mergedBed="sv_calls.bed",
        mergedVCF="sv_calls.vcf",
        mergedVCFgz=expand("{sample}.sv_calls.vcf.gz", sample=config["sample"]),
        mergedBB=expand("{sample}.{op}.bb",sample=config["sample"],op=shortOps),
        inv=expand("{dir}/inversions/inversions.{hap}.bed", dir=config["alnDir"], hap=haps),
        invMerged=expand("inversions.merged.{hap}.bed", hap=haps),
        invHaps=expand("{sample}.inversions.UW.bed",sample=config["sample"]),
        paCov=expand("parent.cov.{pa}.bed", pa=parents),
        isectCounts=expand("{op}.intersect.txt",op=shortOps),
        merged=expand("merged.sources.{op}.txt",op=shortOps),
        bnCompare=["sv_calls.bed.ins.bn.bed", "sv_calls.bed.ins.bn.bed", "sv_calls.bed.by-pb.bn.bed"],
        bnPlot=expand("BioNano.PacBio.Agreement.{sample}.pdf",sample=config["sample"]),
        plot=config["sample"]+".parental_coverage.pdf",
        inh="inheritance.bed",
        hethom="het_hom_summary.txt"
        

rule MakeMSSMSort:
    input:
        mssm="{mssm}"
    output:
        mssmSorted="{mssm}.sorted"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00"
    shell:
        "module unload anaconda; module load python/2.7.3; module load bedtools/latest; bedtools sort -i {input.mssm} | awk '{{ if ($3-$2 >= 50) print;}}' > {output.mssmSorted}"
        
rule MakeMSSMSupport:
    input:
        mssm="{mssm}.sorted",
        reads=config["reads"]
    output:
        cov="{mssm}.cov"
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=04:00:00",
        ref=config["ref"]
    shell:
        "module unload anaconda; module load python/2.7.3; "+SNAKEMAKE_DIR+"/../sv/utils/SpliceVariantsAndCoverageValidate.py --mssm --gaps {input.mssm} --reads {input.reads} --out {output.cov} --nproc 12 --ref {params.ref}"

rule MakeMSSMBN:
    input:
        mssm=lambda wildcards: mssmBaseToPath[wildcards.base] +".sorted"
    output:
        mssmBN="{base}.to-bn.bed"
    params:
        sge_opts=config["sge_small"]
    shell:
        "module unload anaconda && module load python/2.7.3 && ~/projects/HGSVG/hgsvg/sv/utils/mssm/MSSMToUW.py {input.mssm} --keepCoordinates > {output.mssmBN} "
       
   
rule MakeMSSMUW:
    input:
        mssm=lambda wildcards: mssmBaseToPath[wildcards.base] + ".sorted"
    output:
        mssmBN="{base}.to-uw.bed"
    params:
        sge_opts=config["sge_small"]
    shell:
        "module unload anaconda && module load python/2.7.3 && ~/projects/HGSVG/hgsvg/sv/utils/mssm/MSSMToUW.py {input.mssm} > {output.mssmBN} "



rule MakeMSSMBNAnnotation:
    input:
        mssm="{base}.to-uw.bed",
        bionano=config["bionano"]
    output:
        mssmBN="{base}.to-uw.bed.bn"
    params:
        sge_opts=config["sge_small"],
        maxRatio=config["maxRatio"]
    shell:
        "~/projects/HGSVG/hgsvg/sv/utils/mssm/AnnotateBedWithBioNano.py --bionano {input.bionano} --table {input.mssm} --out {output.mssmBN} --source MSSM --maxRatio {params.maxRatio}"

    
rule MakeUWAnnotation:
    input:
        uw=config["uwsv"],
        bionano=config["bionano"]
    output:
        uwBN="uw.bed.bn"
    params:
        sge_opts=config["sge_small"],
        maxRatio=config["maxRatio"]
    shell:
        "~/projects/HGSVG/hgsvg/sv/utils/mssm/AnnotateBedWithBioNano.py --bionano {input.bionano} --table {input.uw} --out {output.uwBN} --source UW --maxRatio {params.maxRatio}"


    
rule FilterCallsByPBAndBNSupportMSSM:
    input:
        mssmCov=lambda wildcards: mssmBaseToPath[wildcards.base] + ".cov",
        mssmToUW="{base}.to-uw.bed",
        mssmToUWBN="{base}.to-uw.bed.bn"
    output:
        filtMSSM="{base}.to-uw.bed.filt"
    params:
        sge_opts=config["sge_small"],
        minPbSupport=config["pbSupport"]
    shell:"""
        echo -e "region\\tnAlt\\tnRef" > {input.mssmToUW}.cov
        cat {input.mssmCov} >> {input.mssmToUW}.cov
        paste {input.mssmToUW} {input.mssmToUW}.cov {input.mssmToUWBN} | sed "s/\\t#/\\t/g" | head -1 > {output.filtMSSM}        
        paste {input.mssmToUW} {input.mssmToUW}.cov {input.mssmToUWBN} | sed "s/\\t#/\\t/g" | bioawk -c hdr '{{ if (((substr($mssmQC,0,4) == "PASS" || $mssmQC == "hap2_resolved_P" || $mssmQC == "hap1_resolved_P") && $nAlt > {params.minPbSupport}) || $bnKey != ".") print;}}' >> {output.filtMSSM}
"""

rule MergeMSSMHaplotypes:
    input:
        mssmh1 = expand("{mssm}.hap1.bed.to-uw.bed.filt",mssm=config["sample"]),
        mssmh2 = expand("{mssm}.hap2.bed.to-uw.bed.filt",mssm=config["sample"]),
        mssmdn = expand("{mssm}.remainingdenovos.bed.to-uw.bed.filt",mssm=config["sample"])
    output:
        mssmMerged = "mssm.bed.filt"
    params:
        sge_opts=config["sge_small"]
    shell:"""
head -1 {input.mssmh1} > {input.mssmh1}.sup
bedtools intersect -v -a {input.mssmdn} -b {input.mssmh1} > {input.mssmdn}.no_h1
cat  {input.mssmh1} {input.mssmdn}.no_h1  | grep -v "^#" | sort -k1,1 -k2,2n -k3,3n >> {input.mssmh1}.sup
bedtools intersect -v -a {input.mssmdn} -b {input.mssmh2} > {input.mssmdn}.no_h2
head -1 {input.mssmh2} > {input.mssmh2}.sup
cat   {input.mssmh2} {input.mssmdn}.no_h2 | grep -v "^#" | sort -k1,1 -k2,2n -k3,3n >> {input.mssmh2}.sup

{SNAKEMAKE_DIR}/../sv/utils/MergeHaplotypes.sh {input.mssmh1}.sup {input.mssmh2}.sup {output.mssmMerged} "svType svLen svSeq qName qStart qEnd mssmQC source region nAlt nRef bnKey bnType bnSize bnRatio bnPB"
"""
    


rule FilterUWCalls:
    input:
        uwsv=config["uwsv"],
        uwbn="uw.bed.bn"
    output:
        uwfilt="uw.bed.filt"
    params:
        sge_opts=config["sge_small"],
        pbSupport=config["pbSupport"]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.uwsv} {input.uwbn} > {output.uwfilt}
paste {input.uwsv} {input.uwbn} | sed "s/\\t#/\\t/g" | bioawk -c hdr '{{ if ($nAlt > {params.pbSupport} || $bnKey != ".") print;}}' | grep -v "^#" | sort -k1,1 -k2,2n -k3,3n | bedtools groupby -c 4 -o first -full | awk  'BEGIN{{OFS="\t";}} NF{{NF-=1}};1;' >> {output.uwfilt}
"""
            
rule MakeBB:
    input:
        bed="{source}.bed.filt"
    output:
        bb=expand("{{source}}.{op}.bb",op=longOps)
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"]
    shell:"""
module load ucsc
egrep "^#|insertion" {input.bed} | bioawk -c hdr 'BEGIN{{OFS="\\t";}} {{ print $_chrom, $tStart, $tEnd, "{wildcards.source}", 1000, "+", $tStart, $tEnd, "0,0,255";}}' | grep -v "^#" | bedtools sort | {SNAKEMAKE_DIR}/../sv/utils/FixCoordinates.py /dev/stdin /dev/stdout {params.ref}.fai >  {wildcards.source}.insertion.bed
bedToBigBed {wildcards.source}.insertion.bed {params.ref}.fai {wildcards.source}.insertion.bb -type=bed9

egrep "^#|deletion" {input.bed} | bioawk -c hdr 'BEGIN{{OFS="\\t";}} {{ print $_chrom, $tStart, $tEnd, "{wildcards.source}", 1000, "+", $tStart, $tEnd, "255,0,0";}}' | grep -v "^#" | bedtools sort | {SNAKEMAKE_DIR}/../sv/utils/FixCoordinates.py /dev/stdin /dev/stdout {params.ref}.fai > {wildcards.source}.deletion.bed
bedToBigBed {wildcards.source}.deletion.bed {params.ref}.fai {wildcards.source}.deletion.bb -type=bed9
"""

rule SeparateByOperation:
    input:
        uw="uw.bed.filt",
        mssm="mssm.bed.filt",
    output:
        svs=expand("{source}.{op}.bed", source=filtered,op=shortOps)
    params:
        sge_opts=config["sge_small"],
    shell:"""
cat {input.uw} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "insertion") print;}}' > {input.uw}.ins.bed
cat {input.uw} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "deletion") print;}}'  > {input.uw}.del.bed
cat {input.mssm} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "insertion") print;}}' > {input.mssm}.ins.bed
cat {input.mssm} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $svType == "deletion") print;}}'  > {input.mssm}.del.bed
"""



rule IntersectMSSMAndUW:
    input:
        svs=expand("{source}.{{op}}.bed", source=filtered)
    output:
        isectCounts="{op}.intersect.txt",
    params:
        sge_opts=config["sge_small"],
    shell:"""
rm -f {output.isectCounts}
for frac in 0.0000001 0.1 0.5; do
  nMSSM=`bedtools intersect -a {filtered[0]}.{wildcards.op}.bed -b {filtered[1]}.{wildcards.op}.bed -r -f $frac -wa | wc -l`
  nUW=`bedtools intersect -b {filtered[0]}.{wildcards.op}.bed -a {filtered[1]}.{wildcards.op}.bed -r -f $frac -wa | wc -l`
  echo -e "$frac\\t$nMSSM\\t$nUW" >> {output.isectCounts}
done
"""

rule CountSpliceClusters:
    input:
        sv="{sv}"
    output:
        svClusters="{sv}.cluster-count"
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        bams=config["reads"],
    shell:
        "{SNAKEMAKE_DIR}/../sv/utils/SpliceVariantsAndCoverageValidate.py --gaps {input.sv} --ref {params.ref} --reads {params.bams} --window 250 --flank 1000 --count {output.svClusters}"


rule SummarizePacBioCovValidation:
    input:
        clusterCount=expand("{uw}.cluster-count", uw=config["uwsv"]) + expand("{mssm}.to-uw.bed.cluster-count",mssm=mssmBaseToPath.keys()),
        mssmCov=expand("{mssm}.to-uw.bed.cov",mssm=mssmBaseToPath.keys()),
        uw=config["uwsv"]
    output:
        clusterSummary="cluster-summary.txt",
    params:
        sge_opts=config["sge_small"],
        minPbSupport=config["pbSupport"],
        sample=config["sample"],
        uw=config["uwsv"]
    shell:"""
nUW=`wc -l {params.uw} | awk '{{ print $1;}}'`
nUWClust=`wc -l {params.uw}.cluster-count | awk '{{ print $1;}}' | awk '{{ print $1;}}'`
nUWPass=`cat {params.uw} | bioawk -c hdr '{{ if ($nAlt > {params.minPbSupport}) print;}}' | wc -l`
echo -e "UW\t$nUW\t$nUWClust\t$nUWPass" > {output.clusterSummary}

nMSSMH1=`wc -l {params.sample}.hap1.bed.to-uw.bed | awk '{{ print $1;}}'`
nMSSMH1Clust=`wc -l {params.sample}.hap1.bed.to-uw.bed.cluster-count | awk '{{ print $1;}}'`
nMSSMH1Pass=`cat {params.sample}.hap1.bed.to-uw.bed.cov | bioawk -c hdr '{{ if ($nAlt > {params.minPbSupport}) print;}}' | wc -l`    
echo -e "MSSMh1\t$nMSSMH1\t$nMSSMH1Clust\t$nMSSMH1Pass" >> {output.clusterSummary}

nMSSMH2=`wc -l {params.sample}.hap2.bed.to-uw.bed | awk '{{ print $1;}}'`
nMSSMH2Clust=`wc -l {params.sample}.hap2.bed.to-uw.bed.cluster-count | awk '{{ print $1;}}' `
nMSSMH2Pass=`cat {params.sample}.hap2.bed.to-uw.bed.cov | bioawk -c hdr '{{ if ($nAlt > {params.minPbSupport}) print;}}' | wc -l`            
echo -e "MSSMh2\t$nMSSMH2\t$nMSSMH2Clust\t$nMSSMH2Pass" >> {output.clusterSummary}

nMSSMDN=`wc -l {params.sample}.remainingdenovos.bed.to-uw.bed | awk '{{ print $1;}}' `
nMSSMDNClust=`wc -l {params.sample}.remainingdenovos.bed.to-uw.bed.cluster-count | awk '{{ print $1;}}' `
nMSSMDNPass=`cat {params.sample}.remainingdenovos.bed.to-uw.bed.cov | bioawk -c hdr '{{ if ($nAlt > {params.minPbSupport}) print;}}' | wc -l`
echo -e "MSSMDN\t$nMSSMDN\t$nMSSMDNClust\t$nMSSMDNPass" >> {output.clusterSummary}
"""

    
rule AnnotateRecirocalOverlap:
    input:
        uw="uw.bed.filt",
        mssm="mssm.bed.filt",
        svs=expand("{source}.{op}.bed", source=filtered,op=shortOps)        
    output:
        annot=expand("{source}.{op}.overlaps.bed", source=filtered,op=shortOps)
    params:
        sge_opts=config["sge_small"],
    shell:"""

echo "ovp" > {input.uw}.ins.overlaps.bed
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.uw}.ins.bed {input.mssm}.ins.bed --index; \
bedtools intersect -a {input.uw}.ins.bed -b {input.mssm}.ins.bed -loj -f 0.5 -r  ; }} | \
bedtools groupby -header  -c 4 -o first -full | \
bioawk -c hdr '{{  if ($chrom_2 != ".") {{ print "SHARED";}} else {{ print "UNIQUE";}} }}'| tail -n +2 >> {input.uw}.ins.overlaps.bed

echo "ovp" > {input.uw}.del.overlaps.bed
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.uw}.del.bed {input.mssm}.del.bed --index; \
bedtools intersect -a {input.uw}.del.bed -b {input.mssm}.del.bed -loj -f 0.5 -r  ; }} | \
bedtools groupby -header  -c 4 -o first -full | \
bioawk -c hdr '{{  if ($chrom_2 != ".") {{ print "SHARED";}} else {{ print "UNIQUE";}} }}'| tail -n +2 >> {input.uw}.del.overlaps.bed


echo "ovp" > {input.mssm}.ins.overlaps.bed
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.mssm}.ins.bed {input.uw}.ins.bed --index; \
bedtools intersect -a {input.mssm}.ins.bed -b {input.uw}.ins.bed -loj -f 0.5 -r  ; }} | \
bedtools groupby -header  -c 4 -o first -full | \
bioawk -c hdr '{{  if ($chrom_2 != ".") {{ print "SHARED";}} else {{ print "UNIQUE";}} }}'| tail -n +2 >> {input.mssm}.ins.overlaps.bed

echo "ovp" > {input.mssm}.del.overlaps.bed
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.mssm}.del.bed {input.uw}.del.bed --index; \
bedtools intersect -a {input.mssm}.del.bed -b {input.uw}.del.bed -loj -f 0.5 -r  ; }} | \
bedtools groupby -header  -c 4 -o first -full | \
bioawk -c hdr '{{  if ($chrom_2 != ".") {{ print "SHARED";}} else {{ print "UNIQUE";}} }}'| tail -n +2 >> {input.mssm}.del.overlaps.bed

"""

#
# Rules for merging callsets:
#
#  1. Pick UW or MSSM calls that overlap with BioNano call and less than 20% size difference.
#  2. From remaining calls, if a UW or MSSM call has overlap with bionano, pick that call. Exclude calls from the other dataset that overlap with the bionnano call.
#  3. Pick remaining UW calls, annotated with which ones are 50% reciprocal overlap with MSSM.
#  4. Pick remaining MSSM calls that are under 10kbp, and not closer than 2kbp from a UW call.
#

rule MakeBionano:
    input:
        sv="{source}.bed.filt.{op}.bed"
    output:
        svbn="{source}.bed.filt.{op}.bed.bn"
    params:
        sge_opts=config["sge_small"]
    shell:"""
bioawk -c hdr -t '{{ if ($bnKey != ".") print;}}' < {input.sv} > {output.svbn}
"""

rule AnnotateOtherBioNanoOverlap:
    input:
        sv="{source}.bed.filt.{op}.bed",
        bn="{other}.bed.filt.{op}.bed.bn"
    output:
        svAnnot="{source}.bed.filt.{op}.{other}.bed.filt-bn.bed"
    params:
        sge_opts=config["sge_small"]
    shell:"""
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.sv} {input.bn} --index ; \
 bedtools intersect -a {input.sv} -b {input.bn} -loj | bedtools groupby -c 4 -o first -full | awk 'BEGIN{{OFS="\\t";}} NF{{NF-=1}};1'; }}  | bioawk -c hdr '{{ if ($chrom_2 != ".") {{ print "{wildcards.other}_BN"; }} else {{ print ".";}} }}' > {output.svAnnot}
"""

rule AnnotateMSSMDistance:
    input:
        mssm="mssm.bed.filt.{op}.bed",
        uw="uw.bed.filt.{op}.bed"
    output:
        svAnnot="mssm.bed.filt.{op}.bed.uw-dist",
    params:
        sge_opts=config["sge_small"]
    shell:"""
echo "distToUW" > {output.svAnnot}
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.mssm} {input.uw} --index --append distToUW ; \
 bedtools closest -header -a {input.mssm} -b {input.uw} -d | bedtools groupby -c 4 -o first -full; }} | bioawk -c hdr '{{ print $distToUW;}}' | tail -n +3  >> {output.svAnnot}
"""

rule AnnotateUWDistance:
    input:
        mssm="mssm.bed.filt.{op}.bed",
        uw="uw.bed.filt.{op}.bed"
    output:
        svAnnot="uw.bed.filt.{op}.bed.mssm-dist",
    params:
        sge_opts=config["sge_small"]
    shell:"""
echo "distToMSSM" > {output.svAnnot}
{{ {SNAKEMAKE_DIR}/../sv/utils/HeaderMod.py --source {input.uw} {input.mssm} --index --append distToMSSM ; \
 bedtools closest -header -a {input.uw} -b {input.mssm} -d | bedtools groupby -c 4 -o first -full; }} | bioawk -c hdr '{{ print $distToMSSM;}}' | tail -n +3 >> {output.svAnnot}
"""

# Merging rules
#
#  1. Pick UW or MSSM calls that overlap with BioNano call and less
#than 20% size difference.
#  2. From remaining calls, if a UW or MSSM call has overlap with
#bionano, pick that call. Exclude calls from the other dataset that
#overlap with the bionnano call.
#  3. Pick remaining UW calls, annotated with which ones are 50%
#reciprocal overlap with MSSM.
#  4. Pick remaining MSSM calls that are under 10kbp, and not closer
#than 2kbp from a UW call.
#  5. Filter out anything  that overlaps an inversion.

rule MergeCallsets:
    input:
        uw="uw.bed.filt",
        mssm="mssm.bed.filt",
        annot=expand("{filt}.{{op}}.overlaps.bed", filt=filtered),
        svbn=expand("{filt}.{{op}}.bed.bn",filt=filtered),
        uwDist="mssm.bed.filt.{op}.bed.uw-dist",
        mssmDist="uw.bed.filt.{op}.bed.mssm-dist"
    output:
        merged="merged.{op}.bed"
    params:
        sge_opts=config["sge_small"],
        minDist=config["minDist"]
    shell:"""


paste uw.bed.filt.{wildcards.op}.bed uw.bed.filt.{wildcards.op}.mssm.bed.filt-bn.bed | head -1 > uw.bed.filt.{wildcards.op}.bed.no-bionano
paste uw.bed.filt.{wildcards.op}.bed uw.bed.filt.{wildcards.op}.mssm.bed.filt-bn.bed | bioawk -c hdr '{{ if ($mssm_BN == ".") print;}}' >> uw.bed.filt.{wildcards.op}.bed.no-bionano

paste mssm.bed.filt.{wildcards.op}.bed mssm.bed.filt.{wildcards.op}.bed.uw-dist | head -1 > mssm.bed.filt.{wildcards.op}.far.from.uw.bed
paste mssm.bed.filt.{wildcards.op}.bed mssm.bed.filt.{wildcards.op}.bed.uw-dist | bioawk -c hdr '{{ if ($bnKey == "." && $svLen < 10000 && $distToUW > {params.minDist}) print; }}' >> mssm.bed.filt.{wildcards.op}.far.from.uw.bed

{SNAKEMAKE_DIR}/../sv/utils/MergeFiles.py --files  mssm.bed.filt.{wildcards.op}.bed.bn uw.bed.filt.{wildcards.op}.bed.no-bionano mssm.bed.filt.{wildcards.op}.far.from.uw.bed --out {output.merged}
"""


rule CountSourceCalls:
    input:
        uw="uw.bed.filt",
        mssm="mssm.bed.filt",
        annot=expand("{filt}.{{op}}.overlaps.bed", filt=filtered),
        svbn=expand("{filt}.{{op}}.bed.bn",filt=filtered),
        uwDist="mssm.bed.filt.{op}.bed.uw-dist",
        mssmDist="uw.bed.filt.{op}.bed.mssm-dist"
    output:
        merged="merged.sources.{op}.txt"
    params:
        sge_opts=config["sge_small"],
        minDist=config["minDist"]
    shell:"""
mssmBN=`bioawk -c hdr '{{ print $svLen;}}' < mssm.bed.filt.{wildcards.op}.bed.bn | stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`
echo -e "MSSM-BN\\t$mssmBN"   > {output.merged}
uwBN=`paste uw.bed.filt.{wildcards.op}.bed uw.bed.filt.{wildcards.op}.mssm.bed.filt-bn.bed | bioawk -c hdr '{{ if ($mssm_BN == "." && $bnRatio <= 0.1  && $bnRatio != ".") print $svLen}}'  | stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`
echo -e "UW-BN\\t$uwBN" >> {output.merged}
uwBNOnly=`paste uw.bed.filt.{wildcards.op}.bed uw.bed.filt.{wildcards.op}.mssm.bed.filt-bn.bed | bioawk -c hdr '{{ if ($nAlt <= 4 && $mssm_BN == "." && $bnRatio <= 0.1  && $bnRatio != ".") print $svLen}}'  | stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`
echo -e "UW-BNOnly\\t$uwBNOnly" >> {output.merged}
uwRem=`paste uw.bed.filt.{wildcards.op}.bed uw.bed.filt.{wildcards.op}.mssm.bed.filt-bn.bed | bioawk -c hdr '{{ if ($nAlt >= 4 && $mssm_BN == "." && ($bnRatio >= 0.1  || $bnRatio == ".")) print $svLen}}'  | stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`
echo -e "uwRem\\t$uwRem" >> {output.merged}
mssmRem=`paste mssm.bed.filt.{wildcards.op}.bed mssm.bed.filt.{wildcards.op}.bed.uw-dist | bioawk -c hdr '{{ if ($bnKey == "." && $svLen < 10000 && $distToUW > {params.minDist}) print $svLen; }}'  | stats.py | tail -n +2 | tr " " "\\t" | cut -f 1-2`
echo -e "MSSM-REM\\t$mssmRem" >> {output.merged}
"""
    
    
rule CombineMergedToPreInv:
    input:
        merged=expand("merged.{op}.bed", op=shortOps),
    output:
        mergedBed="sv_calls.pre-inv.bed",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
module load numpy/1.11.0    
module load pandas
head -1 merged.ins.bed > {output.mergedBed}
cat {input.merged} | grep -v "^#" | awk '{{ if ($3-$2 >= 50) print; }}' | bedtools sort >> {output.mergedBed}
"""

rule ExcludeEventsOverlappingInversions:
    input:
        preInv="sv_calls.pre-inv.bed",
        invHaps="inversions.bed"
    output:
        svCalls="sv_calls.bed"
    params:
        sge_opts=config["sge_small"],
    shell:
        "bedtools intersect -header -v -a {input.preInv} -b {input.invHaps} > {output.svCalls}"

rule CombineMergedToVCF:
    input:
        svCalls="sv_calls.bed"
    output:
        mergedVCF="sv_calls.vcf",
        mergedVCFgz=config["sample"]+".sv_calls.vcf.gz",
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
module load numpy/1.11.0    
module load pandas

{SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.svCalls} --reference {params.ref} --vcf {output.mergedVCF} --sample {params.sample} --type sv --fields NALT nAlt NREF nRef SVANN svAnn SVREP svRep SVCLASS svClass NTR nTR BN bnKey SOURCE source

bgzip -c {output.mergedVCF} > {output.mergedVCFgz}
tabix {output.mergedVCFgz}

"""

rule CombineMergedToTracks:
    input:
        svCalls="sv_calls.bed"
    output:
        mergedBB=expand("{sample}.{op}.bb",sample=config["sample"],op=shortOps)
    params:
        sge_opts=config["sge_small"],
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
module load numpy/1.11.0    
module load pandas

cat merged.del.bed | bioawk -c hdr 'BEGIN{{OFS="\\t";}} {{ print $_chrom, $tStart, $tEnd, $source, 1000, "+", $tStart, $tEnd, "255,0,0";}}' | grep -v "^#" | bedtools sort | {SNAKEMAKE_DIR}/../sv/utils/FixCoordinates.py /dev/stdin /dev/stdout {params.ref}.fai > {params.sample}.del.bed9
bedToBigBed {params.sample}.del.bed9 {params.ref}.fai {params.sample}.del.bb -type=bed9

cat merged.ins.bed | bioawk -c hdr 'BEGIN{{OFS="\\t";}} {{ print $_chrom, $tStart, $tEnd, $source, 1000, "+", $tStart, $tEnd, "0,0,255";}}' | grep -v "^#" | bedtools sort | {SNAKEMAKE_DIR}/../sv/utils/FixCoordinates.py /dev/stdin /dev/stdout {params.ref}.fai > {params.sample}.ins.bed9

bedToBigBed {params.sample}.ins.bed9 {params.ref}.fai {params.sample}.ins.bb -type=bed9
"""

rule MakeInversions:
    input:
        aln=config["alnDir"] + "/alignments.{hap}.bam"
    output:
        inv=config["alnDir"] + "/inversions/inversions.{hap}.bed"
    params:
        alnDir=config["alnDir"],
        sge_opts="-pe serial 8 -l mfree=3G -l h_rt=24:00:00",
        pbs=config["pbs"],
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
module load samtools/latest
mkdir -p {config.alnDir}
samtools view {input.aln} | {params.pbs}/Inversions/screenInversions /dev/stdin {parmas.ref} {output.inv} -w 10000 -r --noClip -j 12
"""

rule MergeInversions:
    input:
        inv=config["alnDir"] + "/inversions/inversions.{hap}.bed"
    output:
        invMerged="inversions.merged.{hap}.bed"
    params:
        sge_opts=config["sge_small"]
    shell:"""
        echo -e "#chrom\\ttStart\\ttEnd" > {output.invMerged}
        bedtools sort -i {input.inv} | bedtools merge >> {output.invMerged}
"""

rule MergeInversionHaplotypes:
    input:
        invs=expand("inversions.merged.{hap}.bed",hap=haps)
    output:
        invHaps=config["sample"] + ".inversions.UW.bed"
    params:
        sge_opts=config["sge_small"]
    shell:"""
{SNAKEMAKE_DIR}/../sv/utils/MergeHaplotypes.sh {input.invs} {output.invHaps} 
"""



def GetParentBams(w):
    if w == 'fa':
        return config["faBams"]
    else:
        return config["moBams"]
    
rule GenotypeParents:
    input:
        svs="sv_calls.bed",
        bams=lambda wildcards: GetParentBams(wildcards.pa)
    output:
        parent="parent.cov.{pa}.bed"
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=08:00:00",
        ref=config["ref"],
        faBam=config["faBams"],
        moBam=config["moBams"]
    shell:"""
module unload anaconda; module load python/2.7.3;

{SNAKEMAKE_DIR}/../sv/utils/SpliceVariantsAndCoverageValidate.py --gaps {input.svs} --reads {input.bams} --out {output.parent} --nproc 12 --ref {params.ref}
"""

rule MakeInheritanceTable:
    input:
        parents=expand("parent.cov.{pa}.bed",pa=parents),
        svcalls="sv_calls.bed"
    output:
        inh="inheritance.bed"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00"
    shell:"""
bioawk -c hdr '{{ print $nAlt"\\t"$nRef"\\t"$hap;}}' < {input.svcalls} > child.support
paste {input.parents} child.support > {output.inh}
"""

rule RenderGenotypes:
    input:
        inh="inheritance.bed"
    output:
        plot=config["sample"]+".parental_coverage.pdf"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00",
        sample=config["sample"],
        pal=config["pal"]
    shell:
        "module load R/latest; Rscript {SNAKEMAKE_DIR}/../plotting/plot_inheritance.R --inh {input.inh} --sample {params.sample} --pal {params.pal}"

rule SummarizeInheritance:
    input:
        svcalls="sv_calls.bed",
        inh="inheritance.bed"
    output:
        hethom="het_hom_summary.txt"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=01:00:00",
        sample=config["sample"],
    shell:"""
nHet=`bioawk -c hdr '{{ if ($hap=="HAP1" || $hap=="HAP0") print; }}' < {input.svcalls} | wc -l`
nHom=`bioawk -c hdr '{{ if ($hap=="HOM") print; }}' < {input.svcalls} | wc -l`
ratio=`echo "scale=3; $nHet/$nHom" | bc`
nHomTP=`cat {input.inh} | {SNAKEMAKE_DIR}/../sv/utils/IncrementHeaderMatches.py | bioawk -c hdr '{{ if ($hap == "HOM" && $nAlt >= 2 && $nAlt_1 >= 2) print; }}' | wc -l`
nHomShouldBeHet=`cat {input.inh} | {SNAKEMAKE_DIR}/../sv/utils/IncrementHeaderMatches.py | bioawk -c hdr '{{ if ($hap == "HOM" && (($nAlt >= 2 && $nAlt_1 < 2) || ($nAlt < 2 && $nAlt_1 >= 2))) print; }}' | wc -l`
nHomFP=`cat {input.inh} | {SNAKEMAKE_DIR}/../sv/utils/IncrementHeaderMatches.py | bioawk -c hdr '{{ if ($hap == "HOM" && $nAlt < 2 && $nAlt_1 < 2) print; }}' | wc -l`
nHetTP=`cat {input.inh} | {SNAKEMAKE_DIR}/../sv/utils/IncrementHeaderMatches.py | bioawk -c hdr '{{ if ( ($hap == "HAP1" || $hap == "HAP0") && ($nAlt >= 2 || $nAlt_1 >=  2)) print; }}' | wc -l`
nHetFP=`cat {input.inh} | {SNAKEMAKE_DIR}/../sv/utils/IncrementHeaderMatches.py | bioawk -c hdr '{{ if ( ($hap == "HAP1" || $hap == "HAP0") && ($nAlt < 2 && $nAlt_1 < 2))  print; }}' | wc -l`
echo -e "#nHom\\tnHet\\tratio\\thomTP\\thetFPHom\\thomFP\\thetTP\\thetFP" > {output.hethom}
echo -e "$nHom\\t$nHet\\t$ratio\\t$nHomTP\\t$nHomShouldBeHet\\t$nHomFP\\t$nHetTP\\t$nHetFP" >> {output.hethom}
paste sv_calls.bed parent.cov.fa.bed parent.cov.mo.bed | {SNAKEMAKE_DIR}/../sv/utils/IncrementHeaderMatches.py | bioawk -c hdr '{{ if ($nAlt_1 < 2  && $nAlt_2 < 2) print;}}' > sv_calls.unconfirmed.bed
"""
        
    

rule MakeBNCompare:
    input:
        bnvcf=config["bnvcf"],
        svcalls="sv_calls.bed"
    output:
        bnCompare=["sv_calls.bed.ins.bn.bed", "sv_calls.bed.ins.bn.bed", "sv_calls.bed.by-pb.bn.bed"],
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=08:00:00",
    shell:
        "{SNAKEMAKE_DIR}/../sv/utils/IntersectBioNanoWithUW.sh {input.bnvcf} {input.svcalls}"

rule MakeBNPlots:
    input:
        bnCompare=["sv_calls.bed.ins.bn.bed", "sv_calls.bed.ins.bn.bed", "sv_calls.bed.by-pb.bn.bed"]
    output:
        bnPlot=expand("BioNano.PacBio.Agreement.{sample}.pdf",sample=config["sample"])
    params:
        sge_opts="-pe serial 8 -l mfree=1G -l h_rt=08:00:00",
        sample=config["sample"]
    shell:
        "module load R/latest; Rscript {SNAKEMAKE_DIR}/../plotting/plot_bionano.R --pbsv sv_calls.bed.by-pb.bn.bed --bndel sv_calls.bed.del.bn.bed --bnins sv_calls.bed.ins.bn.bed --sample {params.sample}"
        
