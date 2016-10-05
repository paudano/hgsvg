PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
HGSVG=/net/eichler/vol5/home/mchaisso/projects/HGSVG/hgsvg
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta
SAMPLE="SAMPLE"
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
  $(DIR)/diploid/deletion/NONE.bed \
  $(DIR)/diploid/insertion/NONE.bed \
  $(DIR)/diploid/insertions.bb \
  $(DIR)/diploid/deletions.bb \
  $(DIR)/diploid/sv_calls.bed \
  $(DIR)/diploid/deletions.annotated.bed \
  $(DIR)/diploid/insertions.annotated.bed \
  $(DIR)/diploid/sv_calls.vcf \
  $(DIR)/diploid/deletions.hap.bed \
  $(DIR)/diploid/insertions.hap.bed \
  $(DIR)/diploid/deletions.hap.bb \
  $(DIR)/diploid/insertions.hap.bb


$(DIR)/hap0/gaps.bed: $(H0SAM)
	mkdir -p $(DIR)/hap0
	$(PBS)/PrintGaps.py $(REF) $(H0SAM) --minContigLength 200000  --minAlignmentLength 10000 --ignoreHP 3 --minDist 200000 --condense 20 --outFile $@ --fractionMasked

$(DIR)/hap1/gaps.bed: $(H1SAM)
	mkdir -p $(DIR)/hap1
	$(PBS)/PrintGaps.py $(REF) $(H1SAM) --minContigLength 200000  --minAlignmentLength 10000  --ignoreHP 3 --minDist 200000 --condense 20 --outFile $@ --fractionMasked

$(DIR)/hap0/insertions.bed: $(DIR)/hap0/gaps.bed
	cd $(DIR)/hap0 && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/$(H0SAM)

$(DIR)/hap1/insertions.bed: $(DIR)/hap1/gaps.bed
	cd $(DIR)/hap1 && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=$(PWD)/$(H1SAM)

$(DIR)/diploid/deletions.bed: $(DIR)/hap0/insertions.bed $(DIR)/hap1/insertions.bed
	mkdir -p $(DIR)/diploid
	bedtools intersect -a $(DIR)/hap0/deletions.bed -b $(DIR)/hap1/deletions.bed -r -f 0.5 -wao | \
    awk '{ if ($$12 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > $(DIR)/diploid/deletions.0.r.bed
	bedtools intersect -b $(DIR)/hap0/deletions.bed -a $(DIR)/hap1/deletions.bed -r -f 0.5 -wao | \
    awk '{ if ($$12 != ".") print; }' | \
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
    awk '{ if ($$12 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > $(DIR)/diploid/insertions.0.r.bed
	bedtools intersect -b $(DIR)/hap0/insertions.bed -a $(DIR)/hap1/insertions.bed -r -f 0.5 -wao | \
    awk '{ if ($$12 != ".") print; }' | \
    cut -f 1-10 | \
    bedtools groupby -c 4 -o first -full > $(DIR)/diploid/insertions.1.r.bed
	bedtools intersect -v -a $(DIR)/hap0/insertions.bed -b $(DIR)/diploid/insertions.0.r.bed -r -f 1.0 > $(DIR)/diploid/insertions.h0.bed
	bedtools intersect -v -a $(DIR)/hap1/insertions.bed -b $(DIR)/diploid/insertions.1.r.bed -r -f 1.0 > $(DIR)/diploid/insertions.h1.bed
	cp $(DIR)/diploid/insertions.0.r.bed $(DIR)/diploid/insertions.hom.bed
	cat $(DIR)/diploid/insertions.h0.bed  | awk '{ print $$0"\tHAP0";}' >> $@
	cat $(DIR)/diploid/insertions.h1.bed  | awk '{ print $$0"\tHAP1";}' >> $@
	cat $(DIR)/diploid/insertions.hom.bed | awk '{ print $$0"\tHOM";}'  >> $@



$(DIR)/diploid/insertion/NONE.bed: $(DIR)/diploid/insertions.bed
	-mkdir $(DIR)/diploid/insertion
	-mkdir $(DIR)/diploid/insertion/rm
	touch $(DIR)/diploid/aln.sam
	touch $(DIR)/diploid/gaps.bed
	# Fake prerequisites for insertions.bed
	cd $(DIR)/diploid && make -t -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertions.bed
#	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertion/L1.bed
	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam insertion/NONE.bed

$(DIR)/diploid/deletion/NONE.bed: $(DIR)/diploid/deletions.bed
	-mkdir deletion
	-mkdir deletion/rm
	touch $(DIR)/diploid/gaps.bed
	touch $(DIR)/diploid/aln.sam
	# Fake prerequisites for insertions.bed
	cd $(DIR)/diploid && make -t -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletions.bed
#	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletion/L1.bed
	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=gaps.bed ALIGNMENTS=aln.sam deletion/NONE.bed

$(DIR)/diploid/insertions.bb: $(DIR)/diploid/insertions.bed
	cut -f 1-3 $^ | bedtools sort  >  "$^".sorted
	module load ucsc; bedToBigBed "$^".sorted $(REF).fai $@

$(DIR)/diploid/deletions.bb: $(DIR)/diploid/deletions.bed
	cut -f 1-3 $^ | bedtools sort  > "$^".sorted
	module load ucsc; bedToBigBed "$^".sorted $(REF).fai $@

$(DIR)/diploid/deletions.annotated.bed:  $(DIR)/diploid/deletion/L1.bed
	$(HGSVG)/sv/utils/CombineAnnotatedEvents.py $@ $(DIR)/diploid/deletion/

$(DIR)/diploid/insertions.annotated.bed: $(DIR)/diploid/insertion/L1.bed
	$(HGSVG)/sv/utils/CombineAnnotatedEvents.py $@ $(DIR)/diploid/insertion/

$(DIR)/diploid/sv_calls.bed: $(DIR)/diploid/insertions.annotated.bed $(DIR)/diploid/deletions.annotated.bed
	cat  $(DIR)/diploid/insertions.annotated.bed $(DIR)/diploid/deletions.annotated.bed > $@

$(DIR)/diploid/sv_calls.vcf:   $(DIR)/diploid/sv_calls.bed
	module load pandas && $(HGSVG)/sv/utils/variants_bed_to_vcf.py --bed $(DIR)/diploid/sv_calls.bed --ref $(REF) --sample $(SAMPLE) --type sv --vcf $@

$(DIR)/diploid/insertions.hap.bed: $(DIR)/diploid/insertions.annotated.bed
	awk '{ col="0,0,255"; if ($$13 == "HAP0") {col="255,0,0";} if ($$13 == "HAP1") {col="0,255,0";} print $$1"\t"$$2"\t"$$3"\t"$$15"\t1000\t+\t"$$2"\t"$$3"\t"col}' <  $^ | sort -k1,1 -k2,2n > $@

$(DIR)/diploid/deletions.hap.bed: $(DIR)/diploid/deletions.annotated.bed
	awk '{ col="0,0,255"; if ($$13 == "HAP0") {col="255,0,0";} if ($$13 == "HAP1") {col="0,255,0";} print $$1"\t"$$2"\t"$$3"\t"$$15"\t1000\t+\t"$$2"\t"$$3"\t"col}' <  $^ | sort -k1,1 -k2,2n > $@


$(DIR)/diploid/insertions.hap.bb: $(DIR)/diploid/insertions.hap.bed
	module load ucsc && bedToBigBed $^ $(REF).fai $@ -type=bed9

$(DIR)/diploid/deletions.hap.bb: $(DIR)/diploid/deletions.hap.bed
	module load ucsc && bedToBigBed $^ $(REF).fai $@ -type=bed9
