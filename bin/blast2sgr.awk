#!/usr/bin/awk -f

BEGIN{
}

NR==FNR{
	len[$1] = $2
}

NR!=FNR{
	found_gene[$2]=1
	if ($9<$10){
		for(i=$9; i<=$10; i++) { hits_plus[$2,i]++ }
	}
	else{
		for(i=$10; i<=$9; i++) { hits_minus[$2,i]++ }
	}
}

END{
	for (j in found_gene){
		if (len [j] ){
			OUTFILE_PLUS = ( FILENAME "_" j ".plus.sgr" )
			OUTFILE_MINUS = ( FILENAME "_" j ".minus.sgr" )
			for (i = 1; i <= len[j]; i++) {
				if ( hits_plus[j,i] ){ printf "%s\t%i\t%i\n", j, i, hits_plus[j,i] > OUTFILE_PLUS }
				else { printf "%s\t%i\t%i\n", j, i, 0 > OUTFILE_PLUS }
				if ( hits_minus[j,i] ){ printf "%s\t%i\t%i\n", j, i, hits_minus[j,i] > OUTFILE_MINUS }
				else { printf "%s\t%i\t%i\n", j, i, 0 > OUTFILE_MINUS }
			}
		}
	}
}
