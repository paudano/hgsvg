PBS=/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts
BLASR=/net/eichler/vol5/home/mchaisso/software/blasr_2/cpp/alignment/bin/blasr
REF=/net/eichler/vol2/eee_shared/assemblies/GRCh38/GRCh38.fasta
MASKER?=repeatmasker
ALIGNMENTS?=alignments.sam
HGSVC=$(HOME)/projects/HGSVG/hgsvg/
GAPS?=gaps.bed

.PRECIOUS:  $(GAPS)

all: $(GAPS) insertions.bed \
		deletions.bed \
	insertion/insertions.fasta\
	deletion/deletions.fasta\
  insertions.bb \
  deletions.bb

#	deletion/full/deletions.partial_masked.bed.trf\
#	insertion/full/insertions.partial_masked.bed.trf
#	insertion/rm/insertions.fasta.out \
#	deletion/rm/deletions.fasta.out \
#	deletion/deletions.annotated.bed \
#	insertion/insertions.annotated.bed \
#	insertion/L1.bed \
#	deletion/L1.bed \
# insertion/NONE.bed\
#  deletion/NONE.bed\

help:
	@echo "Usage: make -f AnnotationPipeline.mak alignments=assembly.sam [INPUT_FASTA=<assemblies_file_name>] "


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
	if [ ! -e insertion/rm ]; then \
    -mkdir -p insertion/rm; \
  fi
	cd insertion && $(HGSVG)/sv/utils/MaskFasta.sh $(MASKER) insertions.fasta rm

deletion/rm/deletions.fasta.out:  deletion/deletions.fasta
	-mkdir -p deletion/rm
	cd deletion && $(HGSVG)/sv/utils/MaskFasta.sh $(MASKER) deletions.fasta rm


MIN_MASK=0.70

#
# The resulting repeat masked files should have mostly sequence annotated as repeat inside of it.
#

insertion/L1.bed: insertion/insertions.annotated.bed
	$(HGSVC)/sv/utils/PrintUniqueEvents.sh insertion/insertions.full_masked.bed insertion

# Partially masked sequences may be fully tandem repeats, annotate this now.

insertion/full/insertions.partial_masked.bed.trf: insertion/insertions.annotated.bed
	$(PBS)/RunTRF.sh insertion/full/insertions.partial_masked.bed

deletion/full/deletions.partial_masked.bed.trf: deletion/deletions.annotated.bed
	$(PBS)/RunTRF.sh deletion/full/deletions.partial_masked.bed


# 
# For remaining annotations, first detect all sequences that are mostly tandem repeat
#   
$(DIR)/diploid/insertion/NONE.bed: $(DIR)/diploid/insertions.bed
	-mkdir -p $(DIR)/diploid/insertion
	-mkdir -p $(DIR)/diploid/insertion/rm
	touch $(DIR)/diploid/aln.sam
	touch $(DIR)/diploid/$(GAPS)
	# Fake prerequisites for insertions.bed
	cd $(DIR)/diploid && make -t -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=$(GAPS) ALIGNMENTS=aln.sam insertions.bed
	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=$(GAPS) ALIGNMENTS=aln.sam insertion/NONE.bed

$(DIR)/diploid/deletion/NONE.bed: $(DIR)/diploid/deletions.bed
	-mkdir -p deletion
	-mkdir -p deletion/rm
	touch $(DIR)/diploid/$(GAPS)
	touch $(DIR)/diploid/aln.sam
	# Fake prerequisites for insertions.bed
	cd $(DIR)/diploid && make -t -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=$(GAPS) ALIGNMENTS=aln.sam deletions.bed
	cd $(DIR)/diploid && make -f $(HGSVG)/stitching/AnnotationPipeline.mak GAPS=$(GAPS) ALIGNMENTS=aln.sam deletion/NONE.bed

$(DIR)/diploid/deletions.annotated.bed: $(DIR)/diploid/deletion/NONE.bed
	$(HGSVG)/sv/utils/CombineAnnotatedEvents.sh $@ $(DIR)/diploid/deletion/

$(DIR)/diploid/insertions.annotated.bed: $(DIR)/diploid/insertion/NONE.bed
	$(HGSVG)/sv/utils/CombineAnnotatedEvents.sh $@ $(DIR)/diploid/insertion/

deletion/L1.bed: deletion/deletions.annotated.bed
	$(HGSVG)/sv/utils/PrintUniqueEvents.sh deletion/deletions.full_masked.bed deletion


insertions.bb: insertions.bed
	cut -f 1-3 $^ | egrep -v "Un|GL|KI" | $(HGSVC)/stitching/TruncateInsertions.py $(REF).fai |  bedtools sort > insertions.lim.bed
	bedToBigBed insertions.lim.bed $(REF).fai $@

deletions.bb: deletions.bed
	cut -f 1-3 $^ | egrep -v "Un|GL|KI" | $(HGSVC)/stitching/TruncateInsertions.py $(REF).fai |  bedtools sort > deletions.lim.bed
	bedToBigBed deletions.lim.bed $(REF).fai $@



clean:
	rm -rf insertion
	rm -rf deletion
	rm -f insertions.bed
	rm -f deletions.bed
