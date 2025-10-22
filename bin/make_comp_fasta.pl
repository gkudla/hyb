#!/usr/bin/perl

my $barcodes=0;

while (<>){
	chomp;
	my ($Fld1,$Fld2) = split("\t", $_, -1);
	if( $.==1 ){
		if( $_ =~ /^\S*_[A-Z]+\s/ ){
#			print "collapsing file, assuming random barcode information encoded in read names\n";
			$barcodes=1;
		}
		else{
#			print "collapsing file, assuming no random barcodes, or random barcodes not processed yet\n";
			$barcodes=0;
		}
	}
	if ( $barcodes==1 ){
		my $read_id = (split(/\s+/, $Fld1))[0];
		@nm = split(/_/, $read_id, -1);
		$barcode = $nm[-1];		# the barcode is the fragment of read name after the last '_'
		$cnt{$Fld2}++;
		if (!$fnd{$Fld2, $barcode}) {
			$fnd{$Fld2, $barcode}++;
			$cnt_barcodes{$Fld2}++;
		}
	}
	else{
		$cnt{$Fld2}++;
	}
}

my $tmp_nm = "$$.tmp";
open (TMP, ">$tmp_nm") || die "Cannot open temporary file: $tmp_nm\n";

if ( $barcodes==1 ){
	foreach $i (keys %cnt) {
		print TMP "$i\t$cnt{$i}\t$cnt_barcodes{$i}\n";
	}
	$syscall='sort -k2,2nr -k3,3nr '.$tmp_nm.' | awk \'{printf ">%s-%s_%s\n%s\n", NR, $3, $2, $1}\'';
	system ($syscall);
	system( "rm $tmp_nm" );
}

else{
	foreach $i (keys %cnt) {
		print TMP "$i\t$cnt{$i}\n";
	}
	$syscall='sort -k2,2nr '.$tmp_nm.' | awk \'{printf ">%s_%s\n%s\n", NR, $2, $1}\'';
	system ($syscall);
	system( "rm $tmp_nm" );
}


