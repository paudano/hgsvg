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

configfile: "config.json"


SLOP_FOR_SV_SEQUENCE_POSITIONS = 5000


faiFile = open(config['ref']+".fai")
chroms = [l.split()[0].rstrip() for l in faiFile]

SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
cwd=os.getcwd()

shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

haps=["h0","h1"]
shortHaps=["0", "1"]
rule all:
    input:
        asmBed      = expand("{base}.{hap}.bam.bed",base=config['alnBase'], hap=haps),
        contigBed   = expand("overlaps/overlaps.{hap}.{chrom}.ctg0.bed",hap=haps,chrom=chroms),
        asmFasta    = expand("{base}.{hap}.bam.fasta", base=config['alnBase'], hap=haps),
        asmOverlaps = expand("overlaps/overlap.{hap}.{chrom}.txt", hap=haps, chrom=chroms),
        asmGraphs   = expand("overlaps/overlap.{hap}.{chrom}.txt.gml", hap=haps, chrom=chroms),
        asmPaths    = expand("overlaps/overlap.{hap}.{chrom}.txt.path", hap=haps, chrom=chroms),
        asmContigs  = expand("contigs/patched.{hap}.{chrom}.fasta", hap=haps, chrom=chroms),
	chrFasta    = expand("contigs.{hap}.fasta", hap=haps),
	chrFastaFai = expand("contigs.{hap}.fasta.fai", hap=haps),
	chrAln      = expand("contigs.{hap}.fasta.sam", hap=haps),
        chrBed      = expand("contigs.{hap}.fasta.sam.bed", hap=haps),
	chrBed6     = expand("contigs.{hap}.fasta.sam.bed6", hap=haps),
	chrBB       = expand("contigs.{hap}.fasta.sam.bb", hap=haps),
        indels      = expand("stitching_hap_gaps/hap{hap}/indels.orig.bed", hap=shortHaps),
        indelVCF    = expand("stitching_hap_gaps/hap{hap}/indels.orig.vcf",hap=shortHaps),
        normVCF     = expand("stitching_hap_gaps/hap{hap}/indels.norm.vcf",hap=shortHaps),
        normBed     = expand("stitching_hap_gaps/hap{hap}/indels.norm.bed",hap=shortHaps),                
        annotation  = "stitching_hap_gaps/diploid/insertions.bed",
        indelBed    ="stitching_hap_gaps/diploid/indels.bed"
        
rule MakeIndels:
    input:
       asmSam="contigs.h{hap}.fasta.sam",
    output:
       indels="stitching_hap_gaps/hap{hap}/indels.orig.bed"
    params:
       sge_opts="-pe serial 2 -l mfree=4G -l h_rt=2:00:00",
       sample=config["sample"],
       sd=SNAKEMAKE_DIR,
       ref=config["ref"],
    shell:"""
mkdir -p stitching_hap_gaps/hap{wildcards.hap}
{params.sd}/../sv/utils/PrintGaps.py {params.ref} {input.asmSam} --minLength 2 --maxLength 49 --ignoreHP 5 --outFile {output.indels}
"""

rule ConvertIndelBedToVCF:
    input:
        indelBed="stitching_hap_gaps/hap{hap}/indels.orig.bed"
    output:
        indelVCF="stitching_hap_gaps/hap{hap}/indels.orig.vcf"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
        sample=config["sample"]
    shell:"""
module unload python/3.5.2; module load python/2.7.3; module load numpy/1.11.0; module load bedtools/latest; module load pandas && {SNAKEMAKE_DIR}/../sv/utils/variants_bed_to_vcf.py --bed {input.indelBed} --ref {params.ref} --sample {params.sample} --type indel --vcf /dev/stdout | bedtools sort -header > {output.indelVCF}
"""


rule NormIndelVCF:
    input:
        indelVCF="stitching_hap_gaps/hap{hap}/indels.orig.vcf"
    output:
        indelNormVCF="stitching_hap_gaps/hap{hap}/indels.norm.vcf"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
    shell:"""
vt normalize -r {params.ref} -o {output.indelNormVCF} {input.indelVCF}
"""

rule NormIndelVCFToBed:
    input:
        indelNormVCF="stitching_hap_gaps/hap{hap}/indels.norm.vcf"
    output:
        indelNormBed="stitching_hap_gaps/hap{hap}/indels.norm.bed"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        ref=config["ref"],
    shell:"""
module unload python/3.5.2; module load python/2.7.3; module load numpy/1.11.0; module load bedtools/latest; module load pandas ; {SNAKEMAKE_DIR}/../sv/utils/variants_vcf_to_bed.py --vcf {input.indelNormVCF} --out {output.indelNormBed}
"""

rule CombineHapIndels:
    input:
        indelNormBed=expand("stitching_hap_gaps/hap{hap}/indels.norm.bed",hap=shortHaps)
    output:
        indelBed="stitching_hap_gaps/diploid/indels.bed"
    params:
        sge_opts="-pe serial 1 -l mfree=2G -l h_rt=04:00:00 -l disk_free=4G",
        sd=SNAKEMAKE_DIR
    shell:"""
{params.sd}/../sv/utils/MergeHaplotypes.sh {input.indelNormBed} {output.indelBed} "svType svLen svSeq" 0.8
"""
    
rule MakeAnnotation:
    input:
       asmSam=expand("contigs.{hap}.fasta.sam",hap=haps)
    output:
        annotation="stitching_hap_gaps/diploid/insertions.bed"
    params:
        sge_opts="-l mfree=1G -pe serial 4 -l h_rt=24:00:00"
    shell:
        "module load numpy/1.11.0 && module load pandas && make -f " + SNAKEMAKE_DIR + "/DiploidAnnotation.mak H0SAM=contigs.h0.fasta.sam H1SAM=contigs.h1.fasta.sam DIR=stitching_hap_gaps -j 2"
    
rule MakeChrAsmFasta:
    input:
        asmContig=expand("contigs/patched.{{hap}}.{chrom}.fasta", chrom=chroms)
    output:
        asmFasta="contigs.{hap}.fasta"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00"
    shell:
        "cat {input.asmContig} > {output.asmFasta}"

rule MakeChrAsmFastaFai:
    input:
        asmFasta="contigs.{hap}.fasta"
    output:
        asmFastaFai="contigs.{hap}.fasta.fai"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00"
    shell:
        "samtools faidx {input.asmFasta}"
    
rule MakeAsmAln:
    input:
        asmContigs=expand("contigs/patched.{{hap}}.{chrom}.fasta.sam", chrom=chroms),
        aln="alignments.{hap}.bam"
    output:
        asmSam="contigs.{hap}.fasta.sam"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00"
    shell:
        "samtools view -H {input.aln} > {output.asmSam}; grep -h -v \"^@\" {input.asmContigs} >> {output.asmSam}"
    
    
rule MakeContigAsmAln:
    input:
        asmFasta="contigs/patched.{hap}.{chrom}.fasta"
    output:
        asmSam=temp("contigs/patched.{hap}.{chrom}.fasta.sam")
    params:
        sge_opts="-l mfree=4G -pe serial 4 -l h_rt=04:00:00",
    	ref=config['ref'],
        sd=SNAKEMAKE_DIR,
        td=TMPDIR
        
    shell:"""
module unload python
module unload python/3.5.2
module load python/2.7.3
module load pysam/0.8.4
module load numpy/1.7.0
module load biopython/latest
    
{params.sd}/MapContigs.py --contigs {input.asmFasta} --ref {params.ref} --tmpdir $TMPDIR --blasr blasr --out {output.asmSam} --nproc 4
"""

rule MakeChrAsmBed:
    input:
        asmSam="contigs.{hap}.fasta.sam"
    output:
        asmBed="contigs.{hap}.fasta.sam.bed"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00",
        hgsvg=SNAKEMAKE_DIR+ "/.."
    shell:
        "samToBed {input.asmSam}  --reportIdentity | bedtools sort > {output.asmBed}"

rule MakeChrAsmBed6:
    input:
        asmBed="contigs.{hap}.fasta.sam.bed"
    output:
        asmBed6="contigs.{hap}.fasta.sam.bed6"
    params:
        sge_opts="-l mfree=1G -pe serial 1 -l h_rt=01:00:00",
	hgsvg=SNAKEMAKE_DIR+ "/.."
    shell:
        "{params.hgsvg}/utils/tracks/SamBedToBed6.py {input.asmBed} {output.asmBed6} "

rule MakeAsmBB:
    input:
        asmBed="contigs.{hap}.fasta.sam.bed6"
    output:
        asmBB="contigs.{hap}.fasta.sam.bb"
    params:
        sge_opts="-l h_rt=01:00:00 -l mfree=1G -pe serial 1",
        ref=config['ref']
    shell:
        "module load ucsc  && bedToBigBed {input.asmBed} {params.ref}.fai {output.asmBB} -type=bed6"


rule MakeAsmContigs:
    input:
        asmPath="overlaps/overlap.{hap}.{chrom}.txt.path",
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt",
        asmFasta="alignments.{hap}.bam.fasta"
    output:
        asmContig="contigs/patched.{hap}.{chrom}.fasta"
    params:
        sge_opts=" -l h_rt=01:00:00 -l mfree=4G -pe serial 1"
    shell:"""
~/projects/HGSVG/hgsvg/stitching/PatchPaths.py {input.asmOverlap} {input.asmFasta} {input.asmPath} {output.asmContig}
"""

rule MakeAsmPaths:
    input:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt",
        asmOverlapGraph="overlaps/overlap.{hap}.{chrom}.txt.gml"
    output:
        asmPath="overlaps/overlap.{hap}.{chrom}.txt.path"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=3G  -pe serial 1"
    shell:"""
module unload python
module unload python/3.5.2        
module load python/2.7.3
module load pysam/0.8.4
module load numpy/1.7.0
module load biopython/latest
module load networkx
~/projects/HGSVG/hgsvg/stitching/OverlapGraphToPaths.py {input.asmOverlap} {input.asmOverlapGraph} {output.asmPath}
"""

        
rule MakeAsmGraphs:
    input:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt"
    output:
        asmOverlapGraph="overlaps/overlap.{hap}.{chrom}.txt.gml"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:"""
module unload python/3.5.2
module load python/2.7.3
module load pysam/0.8.4
module load numpy/1.7.0
module load biopython/latest
module load networkx
~/projects/HGSVG/hgsvg/stitching/OverlapsToGraph.py {input.asmOverlap} --out {output.asmOverlapGraph}
"""


#subworkflow AsmOverlapsWorkflow:
#    snakefile: SNAKEMAKE_DIR +"/MakeAsmOverlaps.Snakefile"
#    workdir:cwd
#

rule SplitAsmOverlaps:
    input:
        bed="overlaps/overlaps.{hap}.{chrom}.ctg0.bed",
    output:
        split=dynamic("overlaps/split_{chrom}/overlaps.{hap}.{chrom}.ctg0.{start}.bed"),
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=1:00:00 -cwd -e /dev/null -o /dev/null  ",
        sd=SNAKEMAKE_DIR,
        nOvp=config['overlapsPerJob']
    shell:"""
mkdir -p overlaps/split_{wildcards.chrom};
module unload python
module load python/2.7.3
module load pysam
module load biopython/latest
module load numpy/1.7.0
        
{params.sd}/SplitBedFile.py --bed {input.bed} --n {params.nOvp} --overlap 5 --base overlaps/split_{wildcards.chrom}/overlaps.{wildcards.hap}.{wildcards.chrom}.ctg0
"""

rule MakeSplitOverlaps:
    input:
        bed="overlaps/split_{chrom}/overlaps.{hap}.{chrom}.ctg0.{start}.bed",
        asm="alignments.{hap}.bam.fasta"
    output:
        splitAsmOverlaps="overlaps/split_{chrom}/overlaps.{hap}.{chrom}.ctg0.{start}.txt"
    params:
        sge_opts="-pe serial 6 -l mfree=2G -l h_rt=6:00:00 -cwd -e /dev/null -o /dev/null ",
        ovps=config["overlapsPerJob"],
        bl=config['blasr'],
        
    shell:"""

module unload python/3.5.2
module load python/2.7.3
module load pysam/0.8.4
module load numpy/1.7.0
module load biopython/latest
echo $PYTHONPATH
which python
mkdir -p $TMPDIR ; mkdir -p overlaps/split_{wildcards.chrom}; /net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg/stitching/OverlapContigsOrderedByBed.py {input.bed} {input.asm} --chrom {wildcards.chrom} --out {output.splitAsmOverlaps} --nproc 12 --tmpdir $TMPDIR --blasr {params.bl}
"""

rule MakeAsmOverlaps:
    input:
        bed="overlaps/overlaps.{hap}.{chrom}.ctg0.bed",
        asm="alignments.{hap}.bam.fasta",
        splitAsmOverlaps=dynamic("overlaps/split_{chrom}/overlaps.{hap}.{chrom}.ctg0.{start}.txt")
    output:
        asmOverlap="overlaps/overlap.{hap}.{chrom}.txt"
    params:
        sge_opts="-pe serial 1 -l mfree=1G -l h_rt=1:00:00 -cwd "
    shell:
        "cat {input.splitAsmOverlaps} > {output.asmOverlap}"

rule MakeContigBed:
    input:
        asmBed = "alignments.{hap}.bam.bed"
    output:
        contigBed = expand("overlaps/overlaps.{{hap}}.{chrom}.ctg0.bed",chrom=chroms)
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        """for c in {} ; do  
        grep \"^$c\t\" {{input.asmBed}} > overlaps/overlaps.{{wildcards.hap}}.$c.bed;
        grep -P "/0\t"  overlaps/overlaps.{{wildcards.hap}}.$c.bed > overlaps/overlaps.{{wildcards.hap}}.$c.ctg0.bed;
        done
        true """.format(" ".join(chroms))
    
rule MakeAsmBed:
    input:
        asmBam = "alignments.{hap}.bam"
    output:
        asmBed = "alignments.{hap}.bam.bed"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        "samtools view {input.asmBam} | samToBed /dev/stdin --reportIdentity > {output.asmBed}"

rule MakeAsmFasta:
    input:
        asmBam = "alignments.{hap}.bam"
    output:
        asmFasta = "alignments.{hap}.bam.fasta"
    params:
        sge_opts="-l h_rt=06:00:00 -l mfree=2G  -pe serial 1"
    shell:
        "samtools view {input.asmBam} | awk '{{ print \">\"$1; print $10;}}' | fold | sed '/^$/d' > {output.asmFasta}; samtools faidx {output.asmFasta}"

