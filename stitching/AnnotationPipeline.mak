PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta
MASKER?=repeatmasker
HGSVC=$(HOME)/projects/HGSVG/hgsvg/
HGSVG=$(HGSVC)
GAPS?=gaps.bed
SAMPLE?=SAMMPLE

.PRECIOUS:  $(GAPS)

all: $(GAPS) insertions.bed \
		deletions.bed \
	insertion/insertions.fasta\
	deletion/deletions.fasta\
  insertions.bb \
  deletions.bb 

help:
	@echo "Usage: make -f AnnotationPipeline.mak GAPS=gaps.bed  "

insertions.bed: $(GAPS) 
	egrep "^#|insertion" $(GAPS) | bedtools sort -header > $@

deletions.bed: $(GAPS)
	egrep "^#|deletion" $(GAPS) | bedtools sort -header > $@


#
# First step to annotation is to create a fasta file.
#

insertion/insertions.fasta: insertions.bed
	-mkdir -p insertion
	$(PBS)/GapBedToFasta.py --seqidx 5 insertions.bed insertion/insertions.fasta

deletion/deletions.fasta: deletions.bed
	-mkdir -p deletion
	$(PBS)/GapBedToFasta.py --seqidx 5 deletions.bed deletion/deletions.fasta

#
# Next, the repeat masked files are generated.
#

insertion/rm/insertions.fasta.out:  insertion/insertions.fasta
	-mkdir -p insertion/rm	
	cd insertion && $(HGSVG)/sv/utils/MaskFasta.sh $(MASKER) insertions.fasta rm

deletion/rm/deletions.fasta.out:  deletion/deletions.fasta
	-mkdir -p deletion/rm
	cd deletion && $(HGSVG)/sv/utils/MaskFasta.sh $(MASKER) deletions.fasta rm


MIN_MASK=0.70
deletion/deletions.annotated.bed: deletion/rm/deletions.fasta.out
	$(PBS)/AnnotateGapBed.py deletions.bed $@ deletion/rm/deletions.fasta.out deletion/rm/deletions.fasta.masked
	cat deletion/deletions.annotated.bed | $(PBS)/PrintFractionMasked.py $(MIN_MASK) lt  > deletion/deletions.partial_masked.bed
	cat deletion/deletions.annotated.bed | $(PBS)/PrintFractionMasked.py $(MIN_MASK) gte > deletion/deletions.full_masked.bed

insertion/insertions.annotated.bed: insertion/rm/insertions.fasta.out
	$(PBS)/AnnotateGapBed.py insertions.bed $@ insertion/rm/insertions.fasta.out insertion/rm/insertions.fasta.masked
	cat insertion/insertions.annotated.bed | $(PBS)/PrintFractionMasked.py $(MIN_MASK) lt > insertion/insertions.partial_masked.bed
	cat insertion/insertions.annotated.bed | $(PBS)/PrintFractionMasked.py $(MIN_MASK) gte > insertion/insertions.full_masked.bed

#
# The resulting repeat masked files should have mostly sequence annotated as repeat inside of it.
#
	$(HGSVC)/sv/utils/PrintUniqueEvents.sh insertion/insertions.full_masked.bed insertion

# Partially masked sequences may be fully tandem repeats, annotate this now.

insertion/insertions.partial_masked.bed.trf: insertion/insertions.annotated.bed
	$(PBS)/RunTRF.sh insertion/insertions.partial_masked.bed

deletion/deletions.partial_masked.bed.trf: deletion/deletions.annotated.bed
	$(PBS)/RunTRF.sh deletion/deletions.partial_masked.bed


# 
# For remaining annotations, first detect all sequences that are mostly tandem repeat
#   
insertion/NONE.bed: insertion/L1.bed insertion/insertions.partial_masked.bed.trf
	cat insertion/insertions.partial_masked.trf.bed |$(HGSVC)/sv/utils/FilterPartialMasked.py 0.70 gte > insertion/TRF.bed
	cat insertion/insertions.partial_masked.trf.bed | $(HGSVC)/sv/utils/FilterPartialMasked.py 0.70 lt >  insertion/not_trf.bed
	egrep "^#|NONE" insertion/not_trf.bed > insertion/NONE.bed
	grep "^#" insertion/not_trf.bed > insertion/partial_annotation.bed
	grep -v NONE insertion/not_trf.bed >> insertion/partial_annotation.bed

insertion/L1.bed: insertion/insertions.annotated.bed
	$(HGSVC)/sv/utils/PrintUniqueEvents.sh insertion/insertions.full_masked.bed insertion

deletion/L1.bed: deletion/deletions.annotated.bed
	$(HGSVC)/sv/utils/PrintUniqueEvents.sh deletion/deletions.full_masked.bed deletion

deletion/NONE.bed: deletion/L1.bed deletion/deletions.partial_masked.bed.trf
	cat deletion/deletions.partial_masked.trf.bed | $(HGSVC)/sv/utils/FilterPartialMasked.py 0.70 gte   > deletion/TRF.bed
	cat deletion/deletions.partial_masked.trf.bed | $(HGSVC)/sv/utils/FilterPartialMasked.py 0.70 lt    > deletion/not_trf.bed

	egrep "^#|NONE" deletion/not_trf.bed > deletion/NONE.bed
	grep "^#" deletion/not_trf.bed > deletion/partial_annotation.bed
	grep -v NONE deletion/not_trf.bed >> deletion/partial_annotation.bed

insertions.bb: insertions.bed
	cut -f 1-3 $^ | egrep -v "Un|GL|KI" | $(HGSVC)/stitching/TruncateInsertions.py $(REF).fai |  bedtools sort > insertions.lim.bed
	bedToBigBed insertions.lim.bed $(REF).fai $@

deletions.bb: deletions.bed
	cut -f 1-3 $^ | egrep -v "Un|GL|KI" | $(HGSVC)/stitching/TruncateInsertions.py $(REF).fai |  bedtools sort > deletions.lim.bed
	bedToBigBed deletions.lim.bed $(REF).fai $@

indels.bed: $(GAPS) 
	$(PBS)/PrintGaps.py $(REF) $(ALIGNMENTS) --outFile $@ --ignoreHP 100 --minLength 2 --maxLength 50


clean:
	rm -rf insertion
	rm -rf deletion
	rm -f insertions.bed
	rm -f deletions.bed
