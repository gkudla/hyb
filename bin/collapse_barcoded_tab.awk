#!/usr/bin/awk -f

# this prints, for each sequence, how many times it was found in the input file, and with how many different barcodes it was found

BEGIN{
	FS="\t"
}

{
        nfields=split($1,nm,"_")
        barcode=nm[nfields]			# the barcode is the fragment of read name after the last '_'
        cnt[$2]++
        if(!fnd[$2,barcode]){
                fnd[$2,barcode]++
                cnt_barcodes[$2]++
        }
}

END{
        for(i in cnt){
                print i, cnt[i], cnt_barcodes[i]
        }
}
