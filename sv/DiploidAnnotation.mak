PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
HGSVG=/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
RE F=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta

	H0SAM?=alignments.h0.sam
H1SAM?=alignments.h1.sam
DIR?=hap_gaps

all: $(DIR)/hap0/gaps.bed \
  $(DIR)/hap1/gaps.bed \
  $(DIR)/hap0/insertions.bed \
  $(DIR)/hap1/insertions.bed \
  $(DIR)/diploid/insertions.bed \
  $(DIR)/diploid/deletions.bed \
  alignments.h0.bam \
  alignments.h1.bam \
  $(DIR)/diploid/deletion/L1.bed \
  $(DIR)/diploid/insertion/L1.bed 

samfiles.fofn: samfiles/chr16.12240000-12300000.sam
	ls samfiles | awk '{ print "samfiles/"$$1;}' > $@

alignments.h0.sam: samfiles.fofn
	$(HGSVG)/sv/FilterSamByHaplotype.py samfiles.fofn

alignments.h0.bam: $(H0SAM)
	mkdir -p /var/tmp/mchaisso
	samtools view -bS $(H1SAM) | samtools sort -T /var/tmp/mchaisso/aln.0.$$PPID -o $@
	samtools index $@

alignments.h1.bam: $(H1SAM)
	mkdir -p /var/tmp/mchaisso
	samtools view -bS $(H1SAM) | samtools sort -T /var/tmp/mchaisso/aln.1.$$PPID -o $@
	samtools index $@

$(DIR)/hap0/gaps.bed: $(H0SAM)
	mkdir -p $(DIR)/hap0
	$(PBS)/PrintGaps.py $(REF) $(H0SAM) --minAlignmentLength 30000 --ignoreHP 3 --minDist 1000 --condense 20 --outFile $@

$(DIR)/hap1/gaps.bed: $(H1SAM)
	mkdir -p $(DIR)/hap1
	$(PBS)/PrintGaps.py $(REF) $(H1SAM) --minAlignmentLength 30000  --ignoreHP 3 --minDist 1000 --condense 20 --outFile $@

$(DIR)/hap0/insertions.bed: $(DIR)/hap0/gaps.bed
	cd $(DIR)/hap0 && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/alignments.h0.bam

$(DIR)/hap1/insertions.bed: $(DIR)/hap1/gaps.bed
	cd $(DIR)/hap1 && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/alignments.h1.bam

$(DIR)/diploid/deletions.bed: $(DIR)/hap0/deletions.bed $(DIR)/hap1/deletions.bed
	mkdir -p $(DIR)/diploid
	bedtools intersect -a $(DIR)/hap0/deletions.bed -b $(DIR)/hap1/deletions.bed -r -f 0.5 -wao | \
    awk '{ if ($$11 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > $(DIR)/diploid/deletions.0.r.bed
	bedtools intersect -b $(DIR)/hap0/deletions.bed -a $(DIR)/hap1/deletions.bed -r -f 0.5 -wao | \
    awk '{ if ($$11 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > $(DIR)/diploid/deletions.1.r.bed
	bedtools intersect -v -a $(DIR)/hap0/deletions.bed -b $(DIR)/diploid/deletions.0.r.bed -r -f 1.0 > $(DIR)/diploid/deletions.h0.bed
	bedtools intersect -v -a $(DIR)/hap1/deletions.bed -b $(DIR)/diploid/deletions.1.r.bed -r -f 1.0 > $(DIR)/diploid/deletions.h1.bed
	cp $(DIR)/diploid/deletions.1.r.bed $(DIR)/diploid/deletions.hom.bed

	cat $(DIR)/diploid/deletions.h0.bed | awk '{ print $$0"\tHAP0";}' > $@
	cat $(DIR)/diploid/deletions.h1.bed | awk '{ print $$0"\tHAP1";}' >> $@
	cat $(DIR)/diploid/deletions.hom.bed | awk '{ print $$0"\tHOM";}' >> $@
#	-rm -f $(DIR)/diploid/deletions.1.r.bed 

$(DIR)/diploid/insertions.bed: $(DIR)/hap0/insertions.bed $(DIR)/hap1/insertions.bed
	mkdir -p $(DIR)/diploid
	bedtools intersect -a $(DIR)/hap0/insertions.bed -b $(DIR)/hap1/insertions.bed -r -f 0.5 -wao | \
    awk '{ if ($$11 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > $(DIR)/diploid/insertions.0.r.bed
	bedtools intersect -b $(DIR)/hap0/insertions.bed -a $(DIR)/hap1/insertions.bed -r -f 0.5 -wao | \
    awk '{ if ($$11 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > $(DIR)/diploid/insertions.1.r.bed
	bedtools intersect -v -a $(DIR)/hap0/insertions.bed -b $(DIR)/diploid/insertions.0.r.bed -r -f 1.0 > $(DIR)/diploid/insertions.h0.bed
	bedtools intersect -v -a $(DIR)/hap1/insertions.bed -b $(DIR)/diploid/insertions.1.r.bed -r -f 1.0 > $(DIR)/diploid/insertions.h1.bed
	cp $(DIR)/diploid/insertions.0.r.bed $(DIR)/diploid/insertions.hom.bed
	cat $(DIR)/diploid/insertions.h0.bed  | awk '{ print $$0"\tHAP0";}' > $@
	cat $(DIR)/diploid/insertions.h1.bed  | awk '{ print $$0"\tHAP1";}' >> $@
	cat $(DIR)/diploid/insertions.hom.bed | awk '{ print $$0"\tHOM";}'  >> $@



$(DIR)/diploid/insertion/L1.bed: $(DIR)/diploid/insertions.bed
	-mkdir $(DIR)/diploid/insertion
	-mkdir $(DIR)/diploid/insertion/rm
	touch $(DIR)/diploid/aln.sam
	touch $(DIR)/diploid/gaps.bed
	cd $(DIR)/diploid && make -t -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertions.bed
	cd $(DIR)/diploid && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertion/L1.bed

$(DIR)/diploid/deletion/L1.bed: $(DIR)/diploid/deletions.bed
	-mkdir deletion
	-mkdir deletion/rm
	touch $(DIR)/diploid/gaps.bed
	touch $(DIR)/diploid/aln.sam
	cd $(DIR)/diploid && make -t -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletions.bed
	cd $(DIR)/diploid && make -f $(PBS)/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletion/L1.bed

$(DIR)/diploid/insertions.hap.bed: $(dir)/diploid/insertions.annotated.bed
	awk '{ col="0,0,255"; if ($$13 == "HAP0") {col="255,0,0";} if ($$13 == "HAP1) {col="0,255,0";} print $$1"\t"$$2\t$$3\t$$15\t1000\t+\t"$$2"\t"$$3"col}' <  $^ > $@

$(DIR)/diploid/deletions.hap.bed: $(dir)/diploid/deletions.annotated.bed
	awk '{ col="0,0,255"; if ($$13 == "HAP0") {col="255,0,0";} if ($$13 == "HAP1) {col="0,255,0";} print $$1"\t"$$2\t$$3\t$$15\t1000\t+\t"$$2"\t"$$3"col}' <  $^ > $@


$(DIR)/diploid/insertions.hap.bb: $(DIR)/diploid/insertions.hap.bed
	module load ucsc && bedToBigBed $^ $(REF).fai $@ -type=bed9

$(DIR)/diploid/deletions.hap.bb: $(DIR)/diploid/deletions.hap.bed
	module load ucsc && bedToBigBed $^ $(REF).fai $@ -type=bed9
