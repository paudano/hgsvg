configfile: "config.json"
import os
import os.path
import subprocess

refFai=open(config['ref']+'.fai')
refChrom = [l.split()[0] for l in refFai]

fofn=open(config['reads'])
inputFiles = [l.rstrip() for l in fofn]
if len(inputFiles) < int(config['chunk']):
    config['chunk'] = str(len(inputFiles))

chunkSize = config['chunk']
    
nChunks = int(len(inputFiles)/chunkSize)
if len(inputFiles) % chunkSize != 0:
    nChunks +=1
    
chunkRange = ["{:02d}".format(i) for i in list(range(1,nChunks))]

alnFiles = [config['aln_base'] + "." + c + ".bam" for c in chunkRange]
SNAKEMAKE_DIR = os.path.dirname(workflow.snakefile)


rule all:
   input:
      fofns=expand("{b}.{r}.fofn", b=config['aln_base'], r=chunkRange),
      bams=expand("bams/{b}.{r}.bam", b=config['aln_base'], r=chunkRange),
      bai=expand("bams/{b}.{r}.bam.bai", b=config['aln_base'], r=chunkRange),
      bamFofn="alignments.txt",
      tagged=expand("tagged/{sample}.{pop}.{chrom}.{aln}.{date}.{dtype}.{phasing}.bam",chrom=refChrom,
                    sample=config['dataset']['sample'],
                    pop=config['dataset']['pop'],
                    aln=config['dataset']['aln'],
                    date=config['dataset']['date'],            
                    dtype=config['dataset']['dtype'],            
                    phasing=config['dataset']['phasing']),
      manifest=expand("manifest/{sample}.{pop}.{chrom}.{aln}.{date}.{dtype}.{phasing}.bam.tsv",chrom=refChrom,
                    sample=config['dataset']['sample'],
                    pop=config['dataset']['pop'],
                    aln=config['dataset']['aln'],
                    date=config['dataset']['date'],            
                    dtype=config['dataset']['dtype'],            
                    phasing=config['dataset']['phasing'])

rule MakeManifest:
    params:
        sample=config['dataset']['sample'],
        pop=config['dataset']['pop'],
        aln=config['dataset']['aln'],
        date=config['dataset']['date'],            
        dtype=config['dataset']['dtype'],            
        phasing=config['dataset']['phasing'],
        fofn=config['reads'],
        part="/net/eichler/vol5/home/mchaisso/projects/pbgreedyphase/partitionByPhasedSNVs",
        vcf=config['vcf'],
        ref=config['ref'],
        sge_opts="-pe serial 1 -l h_rt=24:00:00 -l mfree=1G -N tag"
    input:
        bam=expand("tagged/{sample}.{pop}.{{chrom}}.{aln}.{date}.{dtype}.{phasing}.bam{suffix}",
                    sample=config['dataset']['sample'],
                    pop=config['dataset']['pop'],
                    aln=config['dataset']['aln'],
                    date=config['dataset']['date'],            
                    dtype=config['dataset']['dtype'],            
                    phasing=config['dataset']['phasing'], suffix=["", ".bai"])
    output:
        manifest=expand("manifest/{sample}.{pop}.{{chrom}}.{aln}.{date}.{dtype}.{phasing}.bam{{suffix}}tsv",
                    sample=config['dataset']['sample'],
                    pop=config['dataset']['pop'],
                    aln=config['dataset']['aln'],
                    date=config['dataset']['date'],            
                    dtype=config['dataset']['dtype'],            
                    phasing=config['dataset']['phasing'])
    shell:
        SNAKEMAKE_DIR+"/../utils/uploading/MakeManifest.sh {input.bam} {output.manifest}"
        
      
rule AddTag:
    input:
        bamFofn="alignments.txt"
    params:
        sample=config['dataset']['sample'],
        pop=config['dataset']['pop'],
        aln=config['dataset']['aln'],
        date=config['dataset']['date'],            
        dtype=config['dataset']['dtype'],            
        phasing=config['dataset']['phasing'],
        fofn=config['reads'],
        part="/net/eichler/vol5/home/mchaisso/projects/pbgreedyphase/partitionByPhasedSNVs",
        vcf=config['vcf'],
        ref=config['ref'],
        sge_opts="-pe serial 8 -l h_rt=24:00:00 -l mfree=1G -N tag"
    output:
        bam=expand("tagged/{sample}.{pop}.{{chrom}}.{aln}.{date}.{dtype}.{phasing}.bam",
                    sample=config['dataset']['sample'],
                    pop=config['dataset']['pop'],
                    aln=config['dataset']['aln'],
                    date=config['dataset']['date'],            
                    dtype=config['dataset']['dtype'],            
                    phasing=config['dataset']['phasing'])
    shell:
        "mkdir -p tagged; tabix -h {params.vcf} {wildcards.chrom} > tagged/{wildcards.chrom}.vcf; samtools merge -R {wildcards.chrom} -u -f -b {input.bamFofn} /dev/stdout  | samtools view -h - | {params.part} --sam /dev/stdin --vcf tagged/{wildcards.chrom}.vcf --ref {params.ref} --h1 /dev/stdout --tag TP | samtools view -@ 8 -bS - -o {output.bam}; samtools index {output.bam}"


rule AlignFofn:
   input:
      bams=expand("bams/{b}.{r}.bam", b=config['aln_base'], r=chunkRange)
   output:
      alnFofn="alignments.txt"
   params:
      sge_opts=" -pe serial 1 -l h_rt=1:00:00 -l mfree=1G "
   run:
      outFile = open(output.alnFofn, 'w')
      for bam in input.bams:
          outFile.write(bam + "\n")
      outFile.close()


def GetTemp():
      if "TMPDIR" in os.environ:
          return os.environ["TMPDIR"]
      if 'tmp' in config:
          return config['tmp']


rule IndexBam:
   input:
       bam="bams/{b}.{r}.bam"
   output:
       bai="bams/{b}.{r}.bam.bai"
   params:
       sge_opts="-l h_rt=24:00:00 -l mfree=4G -pe serial 1 -N bai"
   shell:
       "samtools index {input.bam}"
   
rule AlignFile:
   input:
       fofn="{b}.{r}.fofn"
   params:
       ref=config['ref'],
       tmp=GetTemp(),
       sge_opts="-l h_rt=48:00:00 -l mfree=3G -pe serial 12 -N aln.{r}"
   output:
       bam="bams/{b}.{r}.bam"
   shell:
       "mkdir -p {params.tmp}; /net/eichler/vol5/home/mchaisso/projects/blasr/cpp/alignment/bin/blasr {input.fofn} {params.ref} -out /dev/stdout -sam -nproc 12 -insertion 8 -deletion 8 -mismatch 4 -indelRate 3 -minMapQV 20 -advanceExactMatches 10 -maxMatch 25 -sdpTupleSize 13  -sdpMaxAnchorsPerPosition 20 -clipping subread | samtools view -uS - | samtools sort -@8 -m2G  - -T {params.tmp}/aln -o {output} "

rule MakeFofns:
   input:
       fofn=config['reads']
   output:
       expand("{b}.{r}.fofn", b=config['aln_base'], r=chunkRange)
   params:
       base=config['aln_base'],
       linesPerChunk=config['chunk'],
       sge_opts=" -pe serial 1 -l h_rt=1:00:00 -l mfree=1G"
       
   shell:
       "split -l {params.linesPerChunk} -d  {input.fofn} {params.base}.; for f in `ls {params.base}.[0-9]*`; do mv -f $f $f.fofn; done"

       
