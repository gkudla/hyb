#!/usr/bin/perl

# reads .txt/.hyb file or .blast file with hybrids
# .txt/.hyb file format:
# SEQ_ID	SEQ	dG	BIT1_ID	BIT1_st_rd	BIT1_en_rd	BIT1_st_gene	BIT1_en_gene	BIT1_eval	BIT2_ID BIT2_st_rd      BIT2_en_rd      BIT2_st_gene    BIT2_en_gene    BIT2_eval
# outputs sorted table with frequencies of hybrids

$FS = "\t";

print "#$0\n";
foreach $argnum (0 .. $#ARGV) { print "#$ARGV[$argnum]\n"; }

while (<>) {	
	@Fld = split /$FS/,$_;
	if( $ARGV =~ /.txt$|.hyb$/ ) {
		push @arr, $Fld[3], $Fld[9]
	}
	if( $ARGV =~ /.blast$/ && $_ =~ /.+/ ) {
		push @arr, $Fld[1];
	}
	if( $ARGV =~ /.txt$|.hyb$/ || ( $ARGV =~ /.blast$/ && $_ !~ /.+/ ) ) {
		@arr = sort @arr;
		$a = join( "\t", @arr ); 
		$cnt{$a}++;
#		print "$a\n";
		@arr = ();
		$a = "";
	}
}

open SRT, "| sort -nk 3";

foreach (sort keys %cnt){
	if( $_ =~ /\w+/ ){ print SRT $_,"\t",$cnt{$_},"\n"; } 
}

close SRT;
