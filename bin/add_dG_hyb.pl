#!/usr/bin/perl

# example:
# add_dG_hyb.pl file.hyb file.vienna > file_with_dG.hyb

use warnings;
use strict;
my %dGs = ();

my $hyb_nm = $ARGV[0];
my $vienna_nm = $ARGV[1];

local *HYB;
local *VIENNA;

open (VIENNA, "$vienna_nm") || die "Cannot open file: $vienna_nm\n";
while (<VIENNA>){
	chomp;
	if ($_ =~ /^([^_]*_[^_]*)/){
		my $ID=$1;
		$_ = <VIENNA>;
		$_ = <VIENNA>;
		$_ =~ /\t\(([^)]*)\)/;
		$dGs{$ID} = $1;
	}
}
close VIENNA;

open (HYB, "$hyb_nm") || die "Cannot open file: $hyb_nm\n";
LINE: while (<HYB>){
	chomp;
	my @Fld = split("\t", $_, -1);
	if( ($Fld[0] =~ /^#/) || ($Fld[0] !~ /[a-zA-Z0-9]/)){
		my $line = join("\t", @Fld);
		print "$line";
		next LINE;
	}
	if( $dGs{ $Fld[0] } =~ /[0-9]/ ){
		$Fld[2] = $dGs{ $Fld[0] };
	}
	else{
		$Fld[2] = ".";
	}
	my $line = join("\t", @Fld);
	print "$line\n";
}
close HYB;

exit 0;



