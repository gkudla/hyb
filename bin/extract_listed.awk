#!/usr/bin/awk -f

BEGIN{
	print "# extract_selected.awk" "\n" "# listed sequences (file 1) extracted from file 2 together with description"
}

FNR==1 {
	print "# " FILENAME
}

NR==FNR{
	found[$1] = 1
}

NR!=FNR{
	if (found[$1]) {print $0} 
} 
