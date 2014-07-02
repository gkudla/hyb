# @(#)Makefile  2014-07-02  A.J.Travis

#
# Make shared BLAST, BLAT and Bowtie2 databases
#

blast = hOH7.nin
blat = hOH7.2bit hOH7-microRNA.2bit
bowtie2 = hOH7.1.bt2
tab = hOH7.tab
INDEXES = hOH7.nhd  hOH7.nhi  hOH7.nhr  hOH7.nog  hOH7.nsd  hOH7.nsi  hOH7.nsq

TARGETS = $(blast) $(blat) $(bowtie2) $(tab)

all: $(TARGETS)

%.nin: %.fasta
	makeblastdb -in $^ -dbtype nucl -input_type fasta -hash_index -out $* -logfile $@.log

%.2bit: %.fasta
	faToTwoBit $< $@

%.1.bt2: %.fasta
	bowtie2-build $^ $*

%.tab:  %.fasta
	fasta2tab.awk $< > $@

%.fasta: %.fasta.gz
	gunzip -c $< > $@

clean:
	rm -f *.log

clobber: clean
	rm -f $(TARGETS) $(INDEXES)
