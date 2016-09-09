configfile: "phasetag.json"

fai=open(config["ref"]+".fai")
chroms=[l.split()[0] for l in fai ]
bamFofn=open(config["fofn"])
bamList = [l.rstrip() for l in bamFofn]

rule all:
    input:
        expand("tagged/HG00733.PUR.{chrom}.BLASR.20160825.PacBio.10x-phased.bam",chrom=chroms)
    
#        expand("tagged/{sample}.{population}.{chrom}.BLASR.{date}.PacBio.10x-phased.bam", sample=config['sample'], population=config['population'],chrom=chroms, date=config['date'])


rule addtag:
    input:
        bams=bamList
    params:
        sample=config['sample'],
        population=config['population'],
        date=config['date'],
        fofn=config['fofn'],
        part="/net/eichler/vol5/home/mchaisso/projects/pbgreedyphase/partitionByPhasedSNVs",
        vcf=config['vcf'],
        ref=config['ref']
    output:
         "tagged/HG00733.PUR.{chrom}.BLASR.20160825.PacBio.10x-phased.bam"
#        "tagged/{params.sample}.{params.population}.{chrom}.BLASR.{params.date}.PacBio.10x-phased.bam"
    shell:
        "samtools merge {wildcards.chrom} -u -f -b {params.fofn} /dev/stdout  | samtools view -h - | {params.part} --sam /dev/stdin --vcf {params.vcf} --ref {params.ref} --h1 /dev/stdout --tag TP | samtool view -bS - -o {output}; samtools index {output}"
