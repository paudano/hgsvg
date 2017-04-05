SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
configfile: "encode.json"

shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

accFile = open(config["acc"])
accessions = [line.strip() for line in accFile]
paths = {}
for acc in accessions:
    paths[acc] = config["path"]
#
# Read control from separate list 
controlFile = open(config["control"])
controlAccessions = [line.strip() for line in controlFile]


for acc in controlAccessions:
    paths[acc] = config["control_path"]

accessions+=controlAccessions

haps=["0", "1"]
pairs=["1","2"]
rule all:
    input:
        rmdup=expand("hap{hap}/{acc}.rmdup.bam",hap=haps, acc=accessions),
        peaks=expand("hap{hap}/{acc}/{acc}.peaks.narrowPeak", hap=haps, acc=accessions),
        peaksToHg38=expand("hap{hap}/{acc}/{acc}.peaks.narrowPeak.to_hg38.bed", hap=haps, acc=accessions),
        bdg=expand("hap{hap}/{acc}/{acc}_treat_pileup.bdg", hap=haps, acc=accessions),
        bw=expand("hap{hap}/{acc}/{acc}.{hap}_treat_pileup.bw", hap=haps, acc=accessions),
        bams=expand("hap{hap}/{acc}.bam",hap=haps, acc=accessions),
        bai=expand("hap{hap}/{acc}.bam.bai",hap=haps, acc=accessions),        
#        fastq=expand("{path}/{acc}_{pair}.fastq.gz",path=config["path"],acc=accessions,pair=pairs),
        gencode=expand("hap{hap}/gencode.psl",hap=haps)


rule MakeBW:
    input:
        peaks="hap{hap}/{acc}/{acc}_treat_pileup.bdg",
        ref=lambda wildcards: config["hap"+wildcards.hap]
    output:
        bw="hap{hap}/{acc}/{acc}.{hap}_treat_pileup.bw",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=40G -l h_rt=8:00:00"
    shell:"""
module load ucsc
sort -k1,1 -k2,2n -S2G {input.peaks} > {output.bw}.sorted
bedGraphToBigWig {output.bw}.sorted {input.ref}.fai {output.bw}
rm -f {output.bw}.sorted
"""

rule CallPeaksBDG:
    input:
        rmdup="hap{hap}/{acc}.rmdup.bam",
    output:
        peaks="hap{hap}/{acc}/{acc}_treat_pileup.bdg",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=8:00:00",
    shell:"""
module unload numpy
module load numpy/1.8.1
module load setuptools/25.1.1 
module load MACS/2.1.0
mkdir -p hap{wildcards.hap}/{wildcards.acc}
macs2 callpeak -t {input.rmdup} --outdir hap{wildcards.hap}/{wildcards.acc} -n {wildcards.acc} -f BAMPE --bdg -g 2.7e9 
"""


rule CallPeaks:
    input:
        bdg="hap{hap}/{acc}/{acc}_treat_pileup.bdg",
    output:
        peaks="hap{hap}/{acc}/{acc}.peaks.narrowPeak",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=8:00:00",
    shell:"""
module unload numpy
module load numpy/1.8.1
module load setuptools/25.1.1 
module load MACS/2.1.0
mkdir -p hap{wildcards.hap}/{wildcards.acc}
macs2 bdgpeakcall -i {input.bdg}  --outdir hap{wildcards.hap}/{wildcards.acc} -o {wildcards.acc}.peaks.narrowPeak
"""

rule LiftOverPeaks:
    input:
        peaks="hap{hap}/{acc}/{acc}.peaks.narrowPeak",
    output:
        lifted="hap{hap}/{acc}/{acc}.peaks.narrowPeak.to_hg38.bed",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=8G -l h_rt=8:00:00",
        mcutils=config["mcutils"]
    shell:"""
{params.mcutils}/bin/samLiftover ../NA19240.h{wildcards.hap}.to_hg38.sam {input.peaks} {output.lifted}
"""

rule RemoveDuplicates:
    input:
        bam="hap{hap}/{acc}.bam",
        ref=lambda wildcards: config["hap"+wildcards.hap]
    output:
        rmdup="hap{hap}/{acc}.rmdup.bam",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=8:00:00",
    shell:"""
mkdir -p hap{wildcards.hap}
samtools rmdup {input.bam} {output.rmdup} --reference {input.ref}
samtools index {output.rmdup}
"""

rule MapFastq:
    input:
        fastq1=lambda wildcards: paths[wildcards.acc] + wildcards.acc +"_1.fastq.gz",
        fastq2=lambda wildcards: paths[wildcards.acc] + wildcards.acc +"_2.fastq.gz", 
        ref=lambda wildcards: config["hap"+wildcards.hap]
    output:
        bam="hap{hap}/{acc}.bam"
    params:
        sge_opts="-cwd -pe serial 8 -l mfree=3G -l h_rt=8:00:00",
    shell:"""
bwa mem {input.ref} <(zcat {input.fastq1} ) <(zcat {input.fastq2}) -t 12 -w 20 -L 20 |\
samtools view -bS - |\
samtools sort - -T $TMPDIR/{wildcards.acc} -o {output.bam} 
"""

rule MakeBai:
    input:
        bam="hap{hap}/{acc}.bam"
    output:
        bai="hap{hap}/{acc}.bam.bai"    
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=6G -l h_rt=8:00:00"
    shell:
        "samtools index {input.bam}"

rule SplitGencode:
    input:
        gencode=config["gencode"],
    output:
        splitGencode=dynamic("gencode.{index}")
    params:
        nlines=config["blat_lines"],
        sge_opts="-cwd -pe serial 1 -l mfree=1G -l h_rt=1:00:00"
    shell:
        "split --lines {params.nlines} {input.gencode} gencode."


rule MapGencode:
    input:
        gencode="gencode.{index}",
        ref=lambda wildcards: config[wildcards.hap]
    params:
        blat=config["blat"],
        sge_opts="-cwd -pe serial 1 -l mfree=8G -l h_rt=8:00:00"
    output:
        psl=temp("{hap}/gencode.{index}.psl")
    shell:"""
{params.blat} {input.ref} {input.gencode}  {output.psl} -q=dna -t=dna -minScore=200 
"""

rule CombineGencode:
    input:
        pslSplit=dynamic("{hap}/gencode.{index}.psl")
    output:
        psl="{hap}/gencode.psl"
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=1G -l h_rt=8:00:00"
    shell:
        "cat {input.pslSplit} > {output.psl}"

