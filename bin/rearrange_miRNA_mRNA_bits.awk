#!/usr/bin/awk -f

BEGIN{
	OFS="\t"
	GET_ALL_TARGETS_NOT_JUST_mRNA=0
	PRINT_GFF_mRNA=0
	PRINT_GFF_microRNA=0
}

{
        microRNA_first = 0
        microRNA_name = "."
        mRNA_name = "."
        mRNA_start = "."
        mRNA_end = "."
        if (($4 ~ /_microRNA$/ && $10 ~ /_mRNA$/) || ($4 ~ /_microRNA$/ && GET_ALL_TARGETS_NOT_JUST_mRNA) ){
		microRNA_first = 1
		microRNA_name = $4
		microRNA_start = $7
		microRNA_end = $8
		mRNA_name = $10
		mRNA_start = $13
		mRNA_end = $14
	}
        else if ( ($4 ~ /_mRNA$/ && $10 ~ /_microRNA$/) || ($10 ~ /_microRNA$/ && GET_ALL_TARGETS_NOT_JUST_mRNA) ){
		microRNA_first = 0
		microRNA_name = $10
		microRNA_start = $13
		microRNA_end = $14
		mRNA_name = $4
		mRNA_start = $7
		mRNA_end = $8
	}
	else{
		next
	}
}

{
	if( PRINT_GFF_mRNA ){
		if( microRNA_first ){
			bit = 2
		}
		else{
			bit = 1
		}
		printf "%s\t%s\t.\t%i\t%i\t.\t+\t.\t%s\n",mRNA_name,FILENAME,mRNA_start,mRNA_end,$1
#		printf "%s\t%s\t.\t%i\t%i\t.\t+\t.\tgene=%s;seq_ID=%s;bit=%i\n",mRNA_name,FILENAME,mRNA_start,mRNA_end,mRNA_name,$1,bit
	}
	else if( PRINT_GFF_microRNA ){
		if( microRNA_first ){
			bit = 1
		}
		else{
			bit = 2
		}
		printf "%s\t%s\t.\t%i\t%i\t.\t+\t.\t%s\n",microRNA_name,FILENAME,microRNA_start,microRNA_end,$1
#		printf "%s\t%s\t.\t%i\t%i\t.\t+\t.\tgene=%s;seq_ID=%s;bit=%i\n",microRNA_name,FILENAME,microRNA_start,microRNA_end,microRNA_name,$1,bit
	}
	else{
		printf "%s\tmicroRNA_name=%s;mRNA_name=%s;mRNA_start=%i;mRNA_end=%i;microRNA_first=%i\n",$1,microRNA_name,mRNA_name,mRNA_start,mRNA_end,microRNA_first
	}
}

