#!/usr/bin/perl

# Ct2B_GK_3.pl: removed bug that caused the program to crash when a header line started with zero
#		added option HYBRID_SS_MIN to print 5th rather than 8th item from header lines (for hybrid-ss-min)

# Ct2B_GK_2.pl: version modified to always print the nucleotide sequence, even if it's the same as in the previous record

# Ct2B_GK.pl: Modified Zuker's Ct2B.pl program that can handle .ct files produced by hybrid-min (instead of hybrid-ss-min)
# There are two changes: 
# 1. it allows different numbers of nucleotides in each sequence
# 2. it reads sequence names from 8th column in header lines (instead of 5th column)
# G. Kudla, January 13, 2011

# Convert from ct format to Vienna format
# 2 record output - sequence, structure + energy if '=' found in
# header and something follows.
# M. Zuker, February 10, 2009.
#
# INPUT: Use standard input or give file name on the command line.
# OUTPUT: standard output

my $HYBRID_SS_MIN = 0; # set this to 0 to use with hybrid-min output; set to 1 to use with hybrid-ss-min output

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;

use strict;
use warnings;
my (@Bp, @Rec);
my ($energy, $Bp, $Seq);
my $n = 1;
my $curr_line = 1;
my $Last_Seq = '';
while (<>) {
	chomp;
	@Rec = split(' ');
	if ($_ =~ /dG/) {
		$energy = $Seq = $Bp = '';
		$n = $Rec[0];
		my @a = split('=');
		defined($a[1]) && (my @e = split(' ', $a[1]));
		defined($e[0]) && ($energy = "\t(" . $e[0] . ')' );
		print "$Rec[7]\n" if ($HYBRID_SS_MIN==0);     # GK: added this line to print sequence name
		print "$Rec[4]\n" if ($HYBRID_SS_MIN==1);     # GK: added this line to print sequence name
		$curr_line = 1;
	}

	else {
		$Seq .= $Rec[1];
		if ($Rec[4] == 0) {
			$Bp .= '.';
		}
		elsif ($Rec[0] < $Rec[4]) {
			$Bp .= '(';
		}
		else {
			$Bp .= ')';
		}
	}

	if ($curr_line % ($n + 1) == 0) {
			printf "%s\n", $Seq;
			printf "%s%s\n", $Bp,$energy;
	}
#	print "$curr_line\t$_\n"; # for debugging
	$curr_line++;
}
