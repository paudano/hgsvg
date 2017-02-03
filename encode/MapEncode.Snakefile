SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)
configfile: "encode.json"

shell.prefix(". {SNAKEMAKE_DIR}/config.sh; ")

accFile = open(config["acc"])
accessions = [line.strip() for line in accFile]

haps=["0", "1"]
pairs=["1","2"]
rule all:
    input:
        sai=expand("hap{hap}/{acc}_{pair}.sai",hap=haps, acc=accessions, pair=pairs),
        rmdup=expand("hap{hap}/{acc}.rmdup.bam",hap=haps, acc=accessions),
        peaks=expand("hap{hap}/{acc}/{acc}_peaks.narrowPeak", hap=haps, acc=accessions),
        bdg=expand("hap{hap}/{acc}/{acc}_treat_pileup.bdg", hap=haps, acc=accessions),        
        bams=expand("hap{hap}/{acc}.bam",hap=haps, acc=accessions),
        bai=expand("hap{hap}/{acc}.bam.bai",hap=haps, acc=accessions),        
        fastq=expand("{path}/{acc}_{pair}.fastq.gz",path=config["path"],acc=accessions,pair=pairs),
        gencode=expand("hap{hap}/gencode.psl",hap=haps)

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
macs2 callpeak -t {input.rmdup}  --outdir hap{wildcards.hap}/{wildcards.acc} -n {wildcards.acc} -f BAMPE --bdg -g 2.7e9 
"""


rule CallPeaks:
    input:
        bdg="hap{hap}/{acc}/{acc}_treat_pileup.bdg",
    output:
        peaks="hap{hap}/{acc}/{acc}_peaks.narrowPeak",
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=4G -l h_rt=8:00:00",
    shell:"""
module unload numpy
module load numpy/1.8.1
module load setuptools/25.1.1 
module load MACS/2.1.0
mkdir -p hap{wildcards.hap}/{wildcards.acc}
macs2 bdgpeakcall -i {input.bdg}  --outdir hap{wildcards.hap}/{wildcards.acc} -o {wildcards.acc} -f BAMPE --bdg -g 2.7e9 
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
        fastq=config["path"]+"{acc}_{pair}.fastq.gz",
        ref=lambda wildcards: config["hap"+wildcards.hap]
    output:
        sai="hap{hap}/{acc}_{pair}.sai"
    params:
        sge_opts="-cwd -pe serial 8 -l mfree=3G -l h_rt=8:00:00",
    shell:
        "bwa aln {input.ref} <(zcat {input.fastq} ) -f {output.sai}   -n 2 -t 12"

rule MakeBai:
    input:
        bam="hap{hap}/{acc}.bam"
    output:
        bai="hap{hap}/{acc}.bam.bai"    
    params:
        sge_opts="-cwd -pe serial 1 -l mfree=6G -l h_rt=8:00:00"
    shell:
        "samtools index {input.bam}"

rule MakeBam:
    input:
        fastq=[config["path"]+"/{acc}_1.fastq.gz", config["path"]+"/{acc}_2.fastq.gz"],
        sai=["hap{hap}/{acc}_1.sai","hap{hap}/{acc}_2.sai"],
        ref=lambda wildcards: config["hap"+wildcards.hap]
    output:
        bam="hap{hap}/{acc}.bam"
    params:
        sge_opts="-cwd -pe serial 4 -l mfree=6G -l h_rt=8:00:00"
    shell:
        "bwa sampe {input.ref} {input.sai[0]} {input.sai[1]} <(zcat {input.fastq[0]}) <(zcat {input.fastq[1]}) | samtools view -bS - | samtools sort -T $TMPDIR/tmp -o {output.bam}"

rule MapGencode:
    input:
        gencode=config["gencode"],
        ref=lambda wildcards: config[wildcards.hap]
    params:
        blat=config["blat"],
        sge_opts="-cwd -pe serial 8 -l mfree=2G -l h_rt=8:00:00"
    output:
        psl="{hap}/gencode.psl"
    shell:"""
{params.blat} {input.ref} {input.gencode}  {output.psl} -q=dna -t=dna -threads=8 -minScore=200 
"""
    
