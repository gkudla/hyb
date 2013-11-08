#!/usr/bin/awk -f

# inputs tab file with sequences (first) and hybrids txt file (second)
# outputs modified hybrids txt file with sequences in the second column

BEGIN {
	OFS = "\t"
}

NR == FNR{
	seq[$1] = $2
}

NR != FNR{
	if(seq[$1]){
		for ( i = 1; i<=NF; i++ ){
			if ( i==2 ) { printf (seq[$1] "\t") }
			else { printf ($i "\t") } 
		}
		printf "\n"
	}
}

