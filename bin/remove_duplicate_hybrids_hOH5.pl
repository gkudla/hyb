#!/usr/bin/perl

# inputs reference file with microRNAs and mRNAs ranked according to some criteria, and .hyb file
# outputs filtered .hyb file, with exactly one line per seq_ID
# example: filter_hOH5_hyb.pl L1A_L1B_0727A_mtophits.ref data.hyb > filtered.hyb 

# the reference file must have gene names in the first column. other columns are ignored. the order of rows is used for ranking.
# the filtering criteria for selection of hybrids are as follows:
# first, the sum of e-values of the two bits of the hybrid
# second, mRNA-microRNA hybrids are ranked higher than any other hybrids
# third, the rank of the higher-ranked bit in the reference file
# fourth, the rank of lower-ranked bit in the reference file
# finally, the order of lines in the input.hyb file 

use warnings;
use strict;
#use lib '/homes/gkudla/bin';
use Hybrid_long;
my %rank;
my $previous_ID = '';
my @array_hits = ();
my $PREFER_MIM=1;

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;

defined( $ARGV[0] ) && open REF_FILE, "<", $ARGV[0] or die "cannot open file $ARGV[0]";
defined( $ARGV[1] ) && open DATA_FILE, "<", $ARGV[1] or die "cannot open file $ARGV[1]";

LINE:	while( <REF_FILE> ){

	my @Fld = split ( "\t", $_ );
	$rank{ $Fld[0] } = $.;

}

LINE2: while( <DATA_FILE> ){

	my $curr_hit = new Hybrid_long();
	next LINE2 if (($_ =~ /^#/)||!($_ =~ /[a-zA-Z0-9]/));
	$curr_hit -> initialize_hyb( $_, 'default', '.', '.' ) or next LINE2;

	if ( $curr_hit -> seq_ID() ne $previous_ID ){
		print_top_hit( @array_hits ) if (scalar @array_hits) > 0;
		@array_hits = ();
		$previous_ID =  $curr_hit -> seq_ID();
	}

	push (@array_hits, $curr_hit);
}
continue {
        print_top_hit( @array_hits ) if eof;
}
	
sub print_top_hit {
	my $best_hit_so_far = shift @_;
	foreach my $i (@_){
		$best_hit_so_far = rank_two_hybrids( $best_hit_so_far, $i );
	}
	print $best_hit_so_far->line();
	return 0;
#	ATTEMPT AT SORTING THE ARRAY OF HITS, BUT THE CRITERIA WERE TOO COMPLICATED
#	my @sorted_hits = sort { $rank{$a->match_bit_name("_mRNA")} <=> $rank{$b->match_bit_name("_mRNA")} } @array_hits;
#	my $curr_line = ( $sorted_hits[0] -> line() );
#	print $curr_line;
}

sub get_ordered_ranks{	
	my ( $x, $y ) = @_;
	return ( undef, undef ) if ( !defined( $rank{$x} ) && !defined( $rank{$y} ) );
	return ( $rank{ $x }, undef ) if ( !defined( $rank{$y} ) );
	return ( $rank{ $y }, undef ) if ( !defined( $rank{$x} ) );
	return ( $rank{ $x }, $rank{ $y } ) if ( $rank{$x} <= $rank{$y} );
	return ( $rank{ $y }, $rank{ $x } );
}

sub rank_two_hybrids {
# THIS SUBROUTINE RETURNS THE "BETTER" OF TWO HYBRIDS ACCORDING TO e-VALUE, BEING A miRNA-mRNA HYBRID, RANK OF BITS IN THE REFERENCE FILE
# IF THE HYBRIDS DIFFER BY NONE OF THESE CRITERIA, THE SUBROUTINE RETURNS THE FIRST HYBRID
	my ( $x, $y ) = @_;

	my $x_mRNA_name = $x->match_bit_name("_mRNA");
	my $x_miRNA_name = $x->match_bit_name("_microRNA");
	my $y_mRNA_name = $y->match_bit_name("_mRNA");
	my $y_miRNA_name = $y->match_bit_name("_microRNA");

	my ( $x_bit1_nm, $x_bit2_nm ) = $x->get_bit_names();
	my ( $y_bit1_nm, $y_bit2_nm ) = $y->get_bit_names();
	my ( $x_top_rank, $x_bottom_rank ) = get_ordered_ranks( $x_bit1_nm, $x_bit2_nm );
	my ( $y_top_rank, $y_bottom_rank ) = get_ordered_ranks( $y_bit1_nm, $y_bit2_nm );

#	for debugging
#	print "X:\t$x_bit1_nm\t$x_bit2_nm\t" . $rank{$x_bit1_nm} . "\t" . $rank{$x_bit2_nm} . "\t$x_top_rank\t$x_bottom_rank\n";
#	print "Y:\t$y_bit1_nm\t$y_bit2_nm\t" . $rank{$y_bit1_nm} . "\t" . $rank{$y_bit2_nm} . "\t$y_top_rank\t$y_bottom_rank\n";
#	FIRST RANK BY SUM OF E-VALUES
	if ( $x->sum_e_values() != $y->sum_e_values() ){
		$x->sum_e_values() < $y->sum_e_values() ? return $x : return $y;
	}
#	SECOND RANK microRNA-mRNA HYBRIDS HIGHER THAN OTHER HYBRIDS
	if( $PREFER_MIM==1 ){
		return $x if (defined $x_mRNA_name && defined $x_miRNA_name && ((!defined $y_mRNA_name) || (!defined $y_miRNA_name)));
		return $y if (((!defined $x_mRNA_name) || (!defined $x_miRNA_name)) && defined $y_mRNA_name && defined $y_miRNA_name);
	}
#	print "equal on criterion 2\txm:$x_mRNA_name\txmi:$x_miRNA_name\tym:$y_mRNA_name\tymi:$y_miRNA_name\n";
#	NEW THIRD RANK BY RANK OF HIGHER-RANKING BIT	
	if ( defined $x_top_rank && defined $y_top_rank && $x_top_rank != $y_top_rank ){	
		$x_top_rank < $y_top_rank ? return $x : return $y;
	}
	return $x if (defined $x_top_rank && !defined $y_top_rank);
	return $y if (!defined $x_top_rank && defined $y_top_rank);
#	NEW FOURTH RANK BY RANK OF LOWER-RANKING BIT	
	if ( defined $x_bottom_rank && defined $y_bottom_rank && $x_bottom_rank != $y_bottom_rank ){	
		$x_bottom_rank < $y_bottom_rank ? return $x : return $y;
	}
	return $x if (defined $x_bottom_rank && !defined $y_bottom_rank);
	return $y if (!defined $x_bottom_rank && defined $y_bottom_rank);
#	IF NONE OF THE ABOVE THINGS DIFFER (eg for rRNA-rRNA hybrids), THE SUBROUTINE RETURNS THE FIRST HYBRID
	return $x;
}


