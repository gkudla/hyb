#!/usr/bin/awk -f

BEGIN{
	DASHES="-----------------------------------------------------"
	FS="_"
}

/_/{
	print
	split($8,a,"-")
#	miRNA_text=($5 "_microRNA\t" $7 "\t" a[1]) 	# bug fixed 20131014, GK
	miRNA_text=($5 "_" $6 "\t" $7 "\t" a[1])
	miRNA_len=1+a[1]-$7
	mRNA_text=($12 "_" $13 "\t" $14 "\t" $15)
}

/^[actgnACGTN]*$/{
	mRNA_len = length($0) - miRNA_len
	print
	print (substr($0,1,miRNA_len) substr(DASHES,1,mRNA_len) "\t" miRNA_text)
        print (substr(DASHES,1,miRNA_len) substr($0,1+miRNA_len,mRNA_len) "\t" mRNA_text)
}

$1 ~ /^[\(\)\.]*\t/{
	print $0 "\n"
}

