# @(#)Makefile  2014-11-25  A.J.Travis

#
# Make shared BLAST, BLAT and Bowtie2 databases
#

MAKEBLASTDB ?= $(shell which makeblastdb)
FATOTWOBIT ?= $(shell which faToTwoBit)
BOWTIE2BUILD ?= $(shell which bowtie2-build)
FASTA2TAB ?= $(shell which fasta2tab.awk)
GUNZIP ?= $(shell which gunzip)

blast = hOH7.nin
blat = hOH7.2bit hOH7-microRNA.2bit
bowtie2 = hOH7.1.bt2
tab = hOH7.tab
INDEXES = hOH7.nhd  hOH7.nhi  hOH7.nhr  hOH7.nog  hOH7.nsd  hOH7.nsi  hOH7.nsq

TARGETS = $(blast) $(blat) $(bowtie2) $(tab)

all: $(TARGETS)

%.nin: %.fasta
ifneq ($(MAKEBLASTDB),)
	$(MAKEBLASTDB) -in $^ -dbtype nucl -input_type fasta -hash_index -out $* -logfile $@.log
else
	@echo "install BLAST to create $@"
endif

%.2bit: %.fasta
ifneq ($(FATOTWOBIT),)
	$(FATOTWOBIT) $< $@
else
	@echo "install BLAT to create $@"
endif

%.1.bt2: %.fasta
ifneq ($(BOWTIE2BUILD),)
	$(BOWTIE2BUILD) $^ $*
else
	@echo "install Bowtie2 to create $@"
endif

%.tab:  %.fasta
ifneq ($(FASTA2TAB),)
	$(FASTA2TAB) $< > $@
else
	@echo "check hyb/bin is on PATH to create $@"
endif

%.fasta: %.fasta.gz
ifneq ($(GUNZIP),)
	 $(GUNZIP) -c $< > $@
else
	@echo "install gzip to create $@"
endif

# keep intermediate files
.SECONDARY:

clean:
	rm -f *.log

clobber: clean
	rm -f $(TARGETS) $(INDEXES) *.bt2
