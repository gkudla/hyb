#!/usr/bin/perl

my $total_cnt_decompressed = 0;
my %cnt = ();
my %seen = ();

open ( DATA, ("<$ARGV[0]") ) or die "$0: cannot open file $ARGV[0] for reading";
while(<DATA>){
	chomp;
	my @Fld = split ( "\t", $_ );	
	my @nm = split(/_/, $Fld[0]);
	if (!($seen{$Fld[0]}++)){
		$total_cnt_decompressed += $nm[-1];	
	}
}
close DATA;

open ( DATA, ("<$ARGV[0]") ) or die "$0: cannot open file $ARGV[0] for reading";
while (<DATA>){
	chomp;
	my @Fld = split ( "\t", $_ );
	my $gene_nm = $Fld[1];
	my @read_nm = split(/_/, $Fld[0]);
	my $curr_cnt = $read_nm[-1];	
	if( $Fld[9]>$Fld[8] && $_ !~ /_pr-tr|pseudo|different|all_hits/ ){
		$cnt{ $Fld[1] } += $curr_cnt;
	}
}
close DATA;

my @keys = sort { $cnt{$b} <=> $cnt{$a} } keys(%cnt);

foreach $i (@keys){
	$percent=100*$cnt{$i}/$total_cnt_decompressed;
	print "$i\t$cnt{$i}\t$percent\n";
}





