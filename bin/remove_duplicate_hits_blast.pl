#!/usr/bin/perl

# inputs reference file with genes ranked according to some criteria, and .blast file
# outputs filtered .blast file, with exactly one line per seq_ID
# example: remove_duplicate_hits_blast.pl A_E_CACAGC_compressed_hOH6.ref data.blast > filtered.blast 

# the reference file must have gene names in the first column. other columns are ignored. the order of rows is used for ranking.
# the filtering criteria for selection of blast hits are as follows:
# first, the e-values of the hits
# second, miRNA hits are preferred over any other classes of RNA
# third, the rank of hits in the reference file
# finally, the order of lines in the input file 

use warnings;
use strict;
my %rank;
my $previous_ID = '';
my @array_hits = ();

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;

defined( $ARGV[0] ) && open REF_FILE, "<", $ARGV[0] or die "cannot open file $ARGV[0]";
defined( $ARGV[1] ) && open DATA_FILE, "<", $ARGV[1] or die "cannot open file $ARGV[1]";

LINE:	while( <REF_FILE> ){

	my @Fld = split ( "\t", $_ );
	$rank{ $Fld[0] } = $.;

}

LINE2: while( <DATA_FILE> ){

	my @Fld = split ( "\t", $_ );
	my $curr_hit = $_;
	my $seq_ID = $Fld[0]; 

	if ( $seq_ID ne $previous_ID ){
		print_top_hit( @array_hits ) if (scalar @array_hits) > 0;
		@array_hits = ();
		$previous_ID =  $seq_ID;
	}

	push (@array_hits, $curr_hit);
}
continue {
        print_top_hit( @array_hits ) if eof;
}
	
sub print_top_hit {
	my $best_hit_so_far = shift @_;
	foreach my $i (@_){
		$best_hit_so_far = rank_two_hits( $best_hit_so_far, $i );
	}
	print $best_hit_so_far;
	return 0;
}

sub rank_two_hits {
# THIS SUBROUTINE RETURNS THE "BETTER" OF TWO HITS ACCORDING TO e-VALUE, miRNA STATUS, OR RANK OF HITS IN THE REFERENCE FILE
# IF THE HITS DIFFER BY NONE OF THESE CRITERIA, THE SUBROUTINE RETURNS THE FIRST HIT
	my ( $x, $y ) = @_;

	my @X = split ( "\t", $x );
	my @Y = split ( "\t", $y );
	my ($x_gene, $x_evalue) = ($X[1], $X[10]);
	my ($y_gene, $y_evalue) = ($Y[1], $Y[10]);

#	for debugging
#	print "X:\t$x_bit1_nm\t$x_bit2_nm\t" . $rank{$x_bit1_nm} . "\t" . $rank{$x_bit2_nm} . "\t$x_top_rank\t$x_bottom_rank\n";
#	print "Y:\t$y_bit1_nm\t$y_bit2_nm\t" . $rank{$y_bit1_nm} . "\t" . $rank{$y_bit2_nm} . "\t$y_top_rank\t$y_bottom_rank\n";
#	FIRST RANK BY SUM OF E-VALUES
	if ( $x_evalue != $y_evalue ){
		$x_evalue < $y_evalue ? return $x : return $y;
	}
#	SECOND RANK microRNAs HIGHER THAN NON_microRNAs	
	return $x if ($x_gene =~ /_microRNA$/ && $y_gene !~ /_microRNA$/);
	return $y if ($x_gene !~ /_microRNA$/ && $y_gene =~ /_microRNA$/);
#	THIRD RANK BY RANK OF HIGHER-RANKING BIT	
	if ( defined $rank{$x_gene} && defined $rank{$y_gene} && $rank{$x_gene} != $rank{$y_gene} ){	
		$rank{$x_gene} < $rank{$y_gene} ? return $x : return $y;
	}
	return $x if (defined $rank{$x_gene} && !defined $rank{$y_gene});
	return $y if (!defined $rank{$x_gene} && defined $rank{$y_gene});
#	IF NONE OF THE ABOVE THINGS DIFFER, THE SUBROUTINE RETURNS THE FIRST HIT
	return $x;
}


