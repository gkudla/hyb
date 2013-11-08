#!/usr/bin/awk -f

# inputs .tab file with reference sequences and .hyb file (order of input files important)
# outputs two fasta files: XXX.bit_1.fasta and XXX.bit_2.fasta

BEGIN{
	OFS="\t"
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
	miRNA_name = $4
	miRNA_start = $7
	miRNA_end = $8
	mRNA_name = $10
	mRNA_start = $13
	mRNA_end = $14
	# throw away chimeras with antisense bits
	if ( miRNA_end < miRNA_start || mRNA_end < mRNA_start ){
		next
	}
	miRNA_seq = substr( seq[miRNA_name], miRNA_start, 1 + miRNA_end - miRNA_start )
	mRNA_seq = substr( seq[mRNA_name], mRNA_start, 1 + mRNA_end - mRNA_start )
	miRNA_name_extended = ($1 "_" miRNA_name "_" miRNA_start "_" miRNA_end )
	mRNA_name_extended = ($1 "_" mRNA_name "_" mRNA_start "_" mRNA_end )
	print (">" miRNA_name_extended "\n" miRNA_seq) > OUT_miRNA
	print (">" mRNA_name_extended "\n" mRNA_seq) > OUT_mRNA
}



