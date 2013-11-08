#!/usr/bin/awk -f

# inputs .tab file with reference sequences and .hyb file (order of input files important)
# outputs two fasta files: XXX.bit_miRNA.fasta and XXX.bit_mRNA.fasta
# if GET_FULL_miRNA==0, extracts sequences corresponding to miRNA bit, otherwise extracts full miRNA sequences
# mRNA sequences are extended in the 3' direction by NT_EXTEND_mRNA nucleotides

BEGIN{
	OFS="\t"
	GET_FULL_miRNA=1
	NT_EXTEND_mRNA=0
	GET_ALL_TARGETS_NOT_JUST_mRNA=0
}

NR==FNR{
	seq[$1] = $2
}

NR!=FNR && FNR==1{
	match( FILENAME, /.*\./ )
	barename = substr ( FILENAME, 1, RLENGTH-1 )
	OUT_miRNA = (barename ".bit_1.fasta")
	OUT_mRNA = (barename ".bit_2.fasta")
}

NR!=FNR{
	if (($4 ~ /_microRNA$/ && $10 ~ /_mRNA$/) || ($4 ~ /_microRNA$/ && GET_ALL_TARGETS_NOT_JUST_mRNA) ){
		miRNA_name = $4
		miRNA_start = $7
		miRNA_end = $8
		mRNA_name = $10
		mRNA_start = $13
		mRNA_end = $14
	}
	else if ( ($4 ~ /_mRNA$/ && $10 ~ /_microRNA$/) || ($10 ~ /_microRNA$/ && GET_ALL_TARGETS_NOT_JUST_mRNA) ){
		miRNA_name = $10
		miRNA_start = $13
		miRNA_end = $14
		mRNA_name = $4
		mRNA_start = $7
		mRNA_end = $8
	}
	else{
		next
	}
	# throw away chimeras with antisense bits
	if ( miRNA_end < miRNA_start || mRNA_end < mRNA_start ){
		next
	}
	if ( GET_FULL_miRNA ){
		miRNA_start = 1
		miRNA_end = length( seq[ miRNA_name ] )
	}
	if ( NT_EXTEND_mRNA ){
		mRNA_end = mRNA_end + NT_EXTEND_mRNA
	}
	miRNA_seq = substr( seq[miRNA_name], miRNA_start, 1 + miRNA_end - miRNA_start )
	mRNA_seq = substr( seq[mRNA_name], mRNA_start, 1 + mRNA_end - mRNA_start )
	miRNA_name_extended = ($1 "_" miRNA_name "_" miRNA_start "_" miRNA_end )
	mRNA_name_extended = ($1 "_" mRNA_name "_" mRNA_start "_" mRNA_end )
	print (">" miRNA_name_extended "\n" miRNA_seq) > OUT_miRNA
	print (">" mRNA_name_extended "\n" mRNA_seq) > OUT_mRNA
}



