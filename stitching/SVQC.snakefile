configfile: "svqc.json"
haps=config["haps"]
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
ops=["insertion", "deletion"]
dirs=["SVQC", "fill-in"]
shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

rule all:
    input:
        stitching=expand("stitching_hap_gaps/hap{hap}/gaps.bed", hap=haps),
        sortedBed=expand("SVQC/hap{hap}/gaps.sorted",hap=haps),
        reCalled=expand("SVQC/hap{hap}/gaps.recalled",hap=haps),        
        reCalledSorted=expand("SVQC/hap{hap}/gaps.recalled.sorted",hap=haps),        
        indelBed=expand("SVQC/hap{hap}/indels.recalled.bed",hap=haps),
        mergedRetained=expand("SVQC/hap{hap}/indels.svqc.bed",hap=haps),
        svCalls=expand("{d}/diploid/sv_calls.bed",d=dirs),
        combinedTable=expand("{d}/hap{hap}/gaps.bed.support",d=dirs, hap=haps),
        mergedBed="merged/sv_calls.bed",
        mergedVCF="merged/sv_calls.vcf",
	sampleVCF="merged/" + config["sample"] + ".sv_calls.vcf.gz",
        retainedIndels=expand("SVQC/hap{hap}/indels.retained.bed", hap=haps),
        mergedRetainedVCF=expand("SVQC/hap{hap}/indels.svqc.raw.vcf",hap=haps),
        normVCF=expand("SVQC/hap{hap}/indels.svqc.norm.vcf",hap=haps),
        normBED=expand("SVQC/hap{hap}/indels.svqc.norm.bed",hap=haps),
        normOpBed=expand("SVQC/hap{hap}/indels.svqc.norm.{op}.bed",hap=haps,op=ops),
        normOpBedSup=expand("SVQC/hap{hap}/indels.svqc.norm.{op}.bed.support",hap=haps,op=ops),
        dipIndels="SVQC/diploid/indels.bed",
        dipOpIndels=expand("SVQC/diploid/indels.{op}.bed",op=ops),
        dipVCF=expand("SVQC/diploid/{sample}.indels.vcf",sample=config["sample"]),
        recalledRgions=expand("SVQC/hap{hap}/regions.recalled.ref", hap=haps),
        recalledSortedRegions=expand("SVQC/hap{hap}/regions.recalled.ref.sorted", hap=haps),
        localIndel=expand("fill-in/hap{hap}/indels.{op}.bed",hap=haps,op=ops),
        fillInAll=expand("fill-in/hap{hap}/gaps.all_regions.bed",hap=haps),
        fillInHap=expand("fill-in/hap{hap}/gaps.bed",hap=haps),
        fillInCov=expand("fill-in/hap{hap}/gaps.bed.cov",hap=haps)
        
rule SortRecalledRegions:
    input:
        recalled="SVQC/hap{hap}/regions.recalled.ref",
    output:
        rsorted="SVQC/hap{hap}/regions.recalled.ref.sorted",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=1G -l h_rt=04:00:00 -l disk_free=4G",
    shell:
        "bedtools sort -i {input.recalled} > {output.rsorted}"


    
rule CollectRetainedIndels:
    input:
        stitchIndels="stitching_hap_gaps/hap{hap}/indels.bed",
        recalledRegions="SVQC/hap{hap}/regions.recalled.ref.sorted"
    output:
        retainedIndels="SVQC/hap{hap}/indels.retained.bed",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=1G -l h_rt=04:00:00 -l disk_free=4G",
    shell:
        "bedtools intersect -v -a {input.stitchIndels} -b {input.recalledRegions} > {output.retainedIndels} "

    
rule SortGaps:
    input:
        stitching="stitching_hap_gaps/hap{hap}/gaps.bed"
    output:
        svqc="SVQC/hap{hap}/gaps.sorted"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=1G -l h_rt=04:00:00 -l disk_free=4G",
    shell:"""
mkdir -p SVQC/hap{wildcards.hap}
bedtools sort -header -i {input.stitching} > {output.svqc}
"""

rule SplitGaps:
    input:
        gaps="SVQC/hap{hap}/gaps.sorted",
        asm="contigs.h{hap}.fasta"
    output:
        splitGaps=dynamic("SVQC/hap{hap}/split/gaps.bed.{id}")
    params:
        n=config["recall_bin"],
        sge_opts="-cwd -pe serial 1 -l mfree=1G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
        sd=SNAKEMAKE_DIR
    shell:"""
module unload anaconda; module load python/2.7.3; mkdir -p SVQC/hap{wildcards.hap}/split;  {params.sd}/RecallRegionsInGapBed.py --asm {input.asm} --ref {params.ref} --gaps {input.gaps} --split {params.n} --splitDir SVQC/hap{wildcards.hap}/split
"""
    
rule RecallGaps:
    input:
        gaps="SVQC/hap{hap}/split/gaps.bed.{id}",
        asm="contigs.h{hap}.fasta"
    output:
        recalled="SVQC/hap{hap}/split/gaps.recalled.{id}",
        refRegions="SVQC/hap{hap}/split/regions.recalled.ref.{id}",
    params:
        sge_opts="-pe serial 6 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
        ngmlr_cutoff=config["ngmlr_cutoff"],
        sd=SNAKEMAKE_DIR
    shell:"""
module unload anaconda;
module load python/2.7.3;
mkdir -p SVQC/hap{wildcards.hap};
mkdir -p SVQC/hap{wildcards.hap}/indels;
{params.sd}/RecallRegionsInGapBed.py --asm {input.asm} --ref {params.ref} --gaps {input.gaps} --out {output.recalled} --nproc 12 --refRegions SVQC/hap{wildcards.hap}/split/regions.recalled.ref.{wildcards.id} --ngmlr {params.ngmlr_cutoff} --indels {output.recalled}.indel.bed --indelDir SVQC/hap{wildcards.hap}/indels

"""

rule MergeRecalledIndels:
    input:
        gapBed="SVQC/hap{hap}/gaps.recalled"
    output:
        indelBed="SVQC/hap{hap}/indels.recalled.bed"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
    shell:"""
cat SVQC/hap{wildcards.hap}/indels/* | head -1 > {output.indelBed}
cat SVQC/hap{wildcards.hap}/indels/* | grep -v "^#" | bedtools sort >> {output.indelBed}
"""

rule MergeRetainedAndRecalledIndels:
    input:
        indelBed="SVQC/hap{hap}/indels.recalled.bed",
        retainedIndels="SVQC/hap{hap}/indels.retained.bed",
    output:
        mergedRetained="SVQC/hap{hap}/indels.svqc.bed"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
    shell:"""
head -1 {input.indelBed} > {output.mergedRetained}
cat {input.indelBed} {input.retainedIndels} | cut -f 1-10 | grep -v "^#" | bedtools sort >> {output.mergedRetained}
"""

rule ConvertIndelBedToVCF:
    input:
        mergedRetained="SVQC/hap{hap}/indels.svqc.bed"
    output:
        mergedRetainedVCF="SVQC/hap{hap}/indels.svqc.raw.vcf"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
module unload python/3.5.2; module load python/2.7.3; module load numpy/1.11.0; module load bedtools/latest; module load pandas && {SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.mergedRetained} --ref {params.ref} --sample {params.sample} --type indel --vcf /dev/stdout | bedtools sort -header > {output.mergedRetainedVCF}
"""

rule NormIndelVCF:
    input:
        rawVCF="SVQC/hap{hap}/indels.svqc.raw.vcf"
    output:
        normVCF="SVQC/hap{hap}/indels.svqc.norm.vcf"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
    shell:"""
vt normalize -r {params.ref} -o {output.normVCF} {input.rawVCF}
"""

rule NormIndelVCFToBed:
    input:
        vcf="SVQC/hap{hap}/indels.svqc.norm.vcf"
    output:
        bed="SVQC/hap{hap}/indels.svqc.norm.bed"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
    shell:"""
module unload python/3.5.2; module load python/2.7.3; module load numpy/1.11.0; module load bedtools/latest; module load pandas ; {SNAKEMAKE_DIR}/../sv/utils/variants_vcf_to_bed.py --vcf {input.vcf} --out {output.bed}
"""
    
rule SplitNormBed:
    input:
        allbed="SVQC/hap{hap}/indels.svqc.norm.bed"
    output:
        opbed="SVQC/hap{hap}/indels.svqc.norm.{op}.bed"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
    shell:"""
egrep "^#|{wildcards.op}" {input.allbed} > {output.opbed}
"""

rule AddSupport:
    input:
        opbed="SVQC/hap{hap}/indels.svqc.norm.{op}.bed",
        locbed="fill-in/hap{hap}/indels.{op}.bed"
    output:
        opsupport="SVQC/hap{hap}/indels.svqc.norm.{op}.bed.support",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
    shell:"""
cat {input.opbed} | \
 bioawk -c hdr '{{ if (NR==1) {{ print "#oChrom\\toStart\\toEnd\\t"$0;}} else {{ print $1"\\t"$2"\\t"$3"\\t"$0}}}}' |  \
 bedtools slop -header -i stdin -g /net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta.fai  -b 200 | \
 bedtools intersect -header -a stdin -b {input.locbed} -loj  |\
 grep -v "	-1	" |\
 awk '{{ if (NR == 1) {{ print $0"\\tqChrom\\tqStart\\tqEnd\\tqop\\tqsvlen"; }} else {{ print $0;}}}}' | \
  {params.sd}/../indels//CountLocalIndelSupport.py | bioawk -c hdr '{{ if (NR == 1 || $locsup > 1) print;}}' > {output.opsupport}
"""

rule MergeSupport:
    input:
       opsupport=expand("SVQC/hap{{hap}}/indels.svqc.norm.{op}.bed.support",op=ops)
    output:
       mergedopsupport="SVQC/hap{hap}/indels.svqc.norm.bed.support",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
    shell:"""
head -1 {input.opsupport[0]} > {output.mergedopsupport}
cat {input.opsupport} | grep -v "^#" | bedtools sort >> {output.mergedopsupport}
"""

rule HapIndelBedToDiploidBedOp:
    input:
        indels=expand("SVQC/hap{hap}/indels.svqc.norm.{{op}}.bed.support",hap=haps)
    output:
        dipOpIndels="SVQC/diploid/indels.{op}.bed"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../indels/MergeHaplotypes.sh {input.indels} {output.dipOpIndels}
"""


rule HapIndelBedToDiploidBed:
    input:
        dipopindels=expand("SVQC/diploid/indels.{op}.bed",op=ops)
    output:
        dipIndels="SVQC/diploid/indels.bed"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR
    shell:"""
head -1 {input.dipopindels[0]} > {output.dipIndels}
cat {input.dipopindels} | grep -v "^#" | bedtools sort >> {output.dipIndels}
"""

rule IndelBedToVCF:
    input:
        bed="SVQC/diploid/indels.bed"
    output:
        vcf=expand("SVQC/diploid/{sample}.indels.vcf",sample=config["sample"])
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR,
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
module unload python; module load python/2.7.3; module load numpy/1.8.1; module load bedtools/latest; module load pandas/0.20.3 
{SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.bed} --ref {params.ref} --sample {params.sample} --type indel --vcf /dev/stdout | bedtools sort -header > {output.vcf}
        bgzip -c {output.vcf} > {output.vcf}.gz
        tabix {output.vcf}.gz
"""
    
    
    
rule SplicedPBSupport:
    input:
        gaps="{dir}/hap{hap}/split/gaps.recalled.{id}"
    params:
        sge_opts="-pe serial 6 -l mfree=1G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
        bams=config["bams"],
        sd=SNAKEMAKE_DIR
    output:
        pbSupport="{dir}/hap{hap}/split/gaps.recalled_support.{id}"
    shell:"""

module unload anaconda; module load python/2.7.3;
bedtools sort -header -i {input.gaps} > {input.gaps}.tmp
mv -f {input.gaps}.tmp {input.gaps}
mkdir -p SVQC/hap{wildcards.hap};
{params.sd}/../sv/utils/SpliceVariantsAndCoverageValidate.py --gaps {input.gaps} --ref {params.ref} --reads {params.bams} --window 250 --flank 1000 --out {output.pbSupport} --nproc 8
"""        

rule AddSplicedPBSupport:
    input:
        sup="SVQC/hap{hap}/split/gaps.recalled_support.{id}",
        gaps="SVQC/hap{hap}/split/gaps.recalled.{id}"        
    params:
        sge_opts="-pe serial 6 -l mfree=1G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
        bams=config["bams"]
    output:
        pbSupport="SVQC/hap{hap}/split/gaps.gaps_recalled_support.{id}"
    shell:
        "paste  {input.gaps} {input.sup} > {output.pbSupport}"

rule MergeRecallGaps:
    input:
        gaps=dynamic("SVQC/hap{hap}/split/gaps.gaps_recalled_support.{id}"),
        refRegions=dynamic("SVQC/hap{hap}/split/regions.recalled.ref.{id}")
    output:
        recalled="SVQC/hap{hap}/gaps.recalled",
        recalledRegions="SVQC/hap{hap}/regions.recalled.ref"
        
    params:
        sge_opts="-pe serial 6 -l mfree=1G -l h_rt=04:00:00 -l disk_free=4G",
    shell:"""
head -1 SVQC/hap{wildcards.hap}/split/gaps.gaps_recalled_support.0 > {output.recalled}
cat {input.gaps} | grep -v "^#" | bedtools sort >> {output.recalled}
cat {input.refRegions} > {output.recalledRegions}
"""



rule SortRecallGaps:
    input:
        recalled="SVQC/hap{hap}/gaps.recalled"
    output:
        recalledSorted="SVQC/hap{hap}/gaps.recalled.sorted"
    params:
        sge_opts="-pe serial 1 -l mfree=4G -l h_rt=01:00:00 -l disk_free=4G"
    shell:
        "module unload anaconda; module load bedtools/latest; bedtools sort -header -i {input.recalled} > {output.recalledSorted}"

rule RealignGapsMakeDotplots:
    input:
        gaps="SVQC/hap{hap}/gaps.sorted",
        asm="contigs.h{hap}.fasta"
    output:
        dotplots="SVQC/hap{hap}/gaps.bed.realigned.dotplots",
    params:
        sge_opts="-pe serial 6 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"]
    shell:
        "module unload anaconda; module load python/2.7.3; mkdir -p SVQC/hap{wildcards.hap}; mkdir -p SVQC/hap{wildcards.hap}/dotplots; mkdir -p SVQC/hap{wildcards.hap}/indels; " + SNAKEMAKE_DIR + "/RealignRegionsInGapBed.py --asm {input.asm} --ref {params.ref} --gaps {input.gaps} --out /dev/null --nproc 12 --dotplot {output.dotplots} --dotplotDir SVQC/hap{wildcards.hap}/dotplots  --refRegions {output.recalled} "


    
rule MakeSVOpBins:
    input:
        indelbed="SVQC/bins.bed"
    output:
        separateBins=expand("SVQC/bins.{op}.bed", op=ops)
    params:
        sge_opts="-cwd -pe serial 6 -l mfree=1G -l h_rt=12:00:00 -l disk_free=4G",
        ref=config["ref"],
        prindel="/net/eichler/vol5/home/mchaisso/projects/mcst/prindel",
        covWalks=SNAKEMAKE_DIR+"/../sv/utils/CovBinsWalks.py"
    shell:"""
{params.covWalks} SVQC/bins.bed --op ins --out SVQC/bins.insertion.bed;
{params.covWalks} SVQC/bins.bed --op del --out SVQC/bins.deletion.bed;
"""

    
rule MakeCombinedTable:
    input:
        gaps="SVQC/hap{hap}/gaps.recalled.sorted"
    output:
        combined="SVQC/hap{hap}/gaps.bed.support",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=04:00:00 -l disk_free=4G"
    shell:"""
cp -f {input.gaps}  {output.combined}
"""
    
rule AnnotateSVGaps:
    input:
         combined=expand("{{dir}}/hap{hap}/gaps.bed.support",hap=haps)
    output:
         calls="{dir}/diploid/sv_calls.bed"
    params:
        sge_opts="-cwd -pe serial 8 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR
    shell:"""
module unload anaconda; module load python/2.7.3; module load bedtools/latest; module load numpy/1.11.0;
make -f {params.sd}/SVQC.DiploidAnnotation.mak GAPS=gaps.bed.support DIR={wildcards.dir}
"""




rule MakeFillIn:
    input:
        contigSam="alignments.h{hap}.sam",
    output:
        fillInHap="fill-in/hap{hap}/gaps.all.bed"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=04:00:00 -l disk_free=4G",
        pbs=config["pbs"],
        ref=config["ref"],
        sd=SNAKEMAKE_DIR
    shell:"""
mkdir -p fill-in/hap{wildcards.hap}
{params.pbs}/PrintGaps.py {params.ref} {input.contigSam} --maxMasked 10 --minAlignmentLength 30000 --minContigLength 30000 --condense 20  |bedtools sort -header | {params.sd}/../sv/utils/rmdup.py > {output.fillInHap}
"""


rule KeepNotCoveredFillIn:
    input:
        fillInHap="fill-in/hap{hap}/gaps.all.bed",
        contigBed="contigs.h{hap}.fasta.sam.bed"
    output:
        fillIn="fill-in/hap{hap}/gaps.all_regions.bed"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=04:00:00 -l disk_free=4G",
    shell:"""
bedtools intersect -header -v -a {input.fillInHap} -b {input.contigBed} > {output.fillIn}
"""

rule RemoveHeterochromatic:
    input:
        fillInAll="fill-in/hap{hap}/gaps.all_regions.bed"
    output:
        fillIn="fill-in/hap{hap}/gaps.bed"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR
    shell:"""
bedtools intersect -header -a {input.fillInAll} -b {params.sd}/../regions/Regions.Called.bed -wa -u > {output.fillIn}
"""

rule FillInIndel:
    input:
        aln="alignments.h{hap}.bam"
    output:
        indel="fill-in/hap{hap}/indels.{op}.bed",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=6G -l h_rt=04:00:00 ",
        pbs=config["pbs"],
        ref=config["ref"],        
    shell:"""
samtools view -hS {input.aln} | {params.pbs}/PrintGaps.py {params.ref} /dev/stdin --minLength 2 --maxLength 50 | grep {wildcards.op} | bedtools sort -header > {output.indel}
"""
        
        
rule FillInSupport:
    input:
        fillIn="fill-in/hap{hap}/gaps.bed"
    output:
        fillInCov="fill-in/hap{hap}/gaps.bed.cov"
    params:
        sge_opts="-cwd -pe serial 8 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
        bams=config["bams"],
        sd=SNAKEMAKE_DIR
    shell:"""

module unload anaconda; module load python/2.7.3;
{params.sd}/../sv/utils/SpliceVariantsAndCoverageValidate.py --gaps {input.fillIn} --ref {params.ref} --reads {params.bams} --window 250 --flank 1000 --out {output.fillInCov} --nproc 8
"""



rule FilteredFillIn:
    input:
        fillInCov="fill-in/hap{hap}/gaps.bed.cov",
        fillIn="fill-in/hap{hap}/gaps.bed"
    output:
        fillInFilt="fill-in/hap{hap}/gaps.bed.support"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR
    shell:"""
paste {input.fillIn} {input.fillInCov} | bioawk -c hdr '{{ if (substr($0,0,1) == "#" || $nAlt > 3) print;}}' > {output.fillInFilt}
"""

#rule FilteredDiploid:
#    input:
#        fillInFilt=expand("fill-in/hap{hap}/gaps.bed.filt",hap=haps)
#    output:
#        fillInDiploid="fill-in/diploid/sv_calls.bed.support"
#    params:
#        sge_opts="-cwd -pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
#        sd=SNAKEMAKE_DIR
#    shell:"""
#{params.sd}/../sv/utils/MergeHaplotypes.sh {input.fillInFilt} {output.fillInDiploid} "svType svLen svSeq qName qStart qEnd region nAlt nRef"
#"""
#    

rule MakeMergedBed:
    input:
        svqcBed="SVQC/diploid/sv_calls.bed",
        fillinBed="fill-in/diploid/sv_calls.bed"
    output:
        mergedBed="merged/sv_calls.bed"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=04:00:00 -l disk_free=4G",
    shell:"""
module unload anaconda
module load python/2.7.3
module load numpy/1.11.0
module load pandas
mkdir -p merged
bedtools intersect -header -v -a {input.fillinBed} -b {input.svqcBed} > merged/fill-in.bed
head -1 merged/fill-in.bed | awk '{{ print $0"\\tsource";}}' > merged/fill-in.bed.src
tail -n +2 merged/fill-in.bed | awk '{{ print $0"\\tlocal";}}' >>  merged/fill-in.bed.src

head -1 {input.svqcBed} | awk '{{ print $0"\\tsource";}}' > merged/svqc.bed.src 
tail -n +2 {input.svqcBed} |  awk '{{ print $0"\\tstitching";}}' >> merged/svqc.bed.src

 {SNAKEMAKE_DIR}/../sv/utils/Select.py --table merged/fill-in.bed.src  --out merged/fill-in.bed.subset --cols `head -1 {input.svqcBed}` source
 {SNAKEMAKE_DIR}/../sv/utils/MergeFiles.py --files merged/svqc.bed.src merged/fill-in.bed.src | bedtools sort -header > {output.mergedBed}
    
"""
#    .format(SNAKEMAKE_DIR, SNAKEMAKE_DIR)


    
rule MakeMergedVCF:
    input:
        mergedBed="merged/sv_calls.bed"
    output:
        mergedVCF="merged/sv_calls.vcf"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"], 
        sample=config["sample"]
    shell:
        "module unload anaconda; module load python/2.7.3; module load numpy/1.11.0; module load bedtools/latest; module load pandas && {SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.mergedBed} --ref {params.ref} --sample {params.sample} --type sv --vcf /dev/stdout --fields NALT nAlt NREF nRef SRC source | bedtools sort -header > {output.mergedVCF}"


rule MakeSampleVCF:
    input:
        mergedVCF="merged/sv_calls.vcf"
    output:
        sampleVCF="merged/"+config["sample"]+".sv_calls.vcf.gz"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=04:00:00 -l disk_free=4G",
    shell:
        "bgzip -c {input.mergedVCF} > {output.sampleVCF}; tabix {output.sampleVCF}"
        
    
