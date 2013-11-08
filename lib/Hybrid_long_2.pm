#class Hybrid_long
package Hybrid_long;
use strict;

#constructor
sub new {
	my ($class) = @_;
	my $self = {
		_seq_ID		=> undef,
 		_seq		=> undef,
 		_dG		=> undef,
 		_bit1_nm	=> undef,
 		_bit1_st_in_rd	=> undef,
 		_bit1_end_in_rd	=> undef,
 		_bit1_st	=> undef,
 		_bit1_end	=> undef,
 		_bit1_eval	=> undef,
 		_bit2_nm	=> undef,
 		_bit2_st_in_rd	=> undef,
 		_bit2_end_in_rd	=> undef,
 		_bit2_st	=> undef,
 		_bit2_end	=> undef,
 		_bit2_eval	=> undef,
 		_count		=> undef,
		_experiment	=> undef,
		_found_overlap	=> undef,
		_twoway_overlap	=> undef,
		_note		=> undef,
		_line		=> undef,
		_sorted_bit_nm	=> undef,
		_seq_ID_list	=> []
	};
	bless $self, $class;
	return $self;
}

#accessor method for Hybrid_long elements
sub initialize_hyb {
	my ( $self, $line, $exp ) = @_;
	my ( $seq_ID, $seq, $dG, $bit1_nm, $bit1_st_in_rd, $bit1_end_in_rd, $bit1_st, $bit1_end, $bit1_eval, $bit2_nm, $bit2_st_in_rd, $bit2_end_in_rd, $bit2_st, $bit2_end, $bit2_eval, $last_column) = split("\t", $line, -1);
	
	$self->{_seq_ID} = $seq_ID if defined( $seq_ID );
	$self->{_seq} = $seq if defined( $seq );
	$self->{_dG} = $dG if defined( $dG );
	$self->{_bit1_nm} = $bit1_nm if defined( $bit1_nm );
	$self->{_bit1_st_in_rd} = $bit1_st_in_rd if defined( $bit1_st_in_rd );
	$self->{_bit1_end_in_rd} = $bit1_end_in_rd if defined( $bit1_end_in_rd );
	$self->{_bit1_st} = $bit1_st if defined( $bit1_st );
	$self->{_bit1_end} = $bit1_end if defined( $bit1_end );
	$self->{_bit1_eval} = $bit1_eval if defined( $bit1_eval );
	$self->{_bit2_nm} = $bit2_nm if defined( $bit2_nm );
	$self->{_bit2_st_in_rd} = $bit2_st_in_rd if defined( $bit2_st_in_rd );
	$self->{_bit2_end_in_rd} = $bit2_end_in_rd if defined( $bit2_end_in_rd );
	$self->{_bit2_st} = $bit2_st if defined( $bit2_st );
	$self->{_bit2_end} = $bit2_end if defined( $bit2_end );
	$self->{_bit2_eval} = $bit2_eval if defined( $bit2_eval );

	$self->{_line} = $line if defined( $line );
	$self->{_experiment} = $exp if defined( $exp );
	$self->{_found_overlap} = 1;
	$self->{_twoway_overlap} = 0;
	@{ $self->{_seq_ID_list} } = ($seq_ID) ;

	if ( $bit1_nm lt $bit2_nm ){
		$self->{_sorted_bit_nm} = "$bit1_nm\t$bit2_nm";
	}
	else{
		$self->{_sorted_bit_nm} = "$bit2_nm\t$bit1_nm";
	}

	if ( defined( $last_column ) && $last_column =~ /^[0-9]+$/ ) {
		$self->{_count} = $last_column;
	}
	elsif ( defined( $last_column ) && $last_column =~ /^[a-zA-Z_0-9]+=.*/ ){
		$self->{_count} = $2 if ( $last_column =~ /(^|;)count_total=([0-9]+)/ );
		$self->{_twoway_overlap} = $2 if ( $last_column =~ /(^|;)two_way_merged=([0-9]+)/ );
		@{ $self->{_seq_ID_list} } = split(";",$2) if ( $last_column =~ /(^|;)seq_IDs_in_cluster=([^;]+)/ );
	}
	else {
		$self->{_count} = 1;
	}
	return 1;
}

sub same_endpoints_same_strand {
# returns 1 if two target sites represent the same exon-exon junction, 0 otherwise (FOR FINDING INTRONS)
	my ($self, $other) = @_;
	if ( $self->{_bit1_nm} eq $other->{_bit1_nm}     &&     $self->{_bit2_nm} eq $other->{_bit2_nm} ){
		if ( ($self->{_bit1_end} == $other->{_bit1_end})  &&  ($self->{_bit2_st} == $other->{_bit2_st}) ){
			return 1;
		}
	}
	return 0;
}

sub overlaps {
# returns 1 if two target sites overlap, 0 otherwise
	my ($self, $other) = @_;
	my $bit1_ovlp = 0;
	my $bit2_ovlp = 0;
	if ( $self->{_bit1_nm} eq $other->{_bit1_nm}     &&     $self->{_bit2_nm} eq $other->{_bit2_nm} ){
		if ( $self->{_bit1_st} < $self->{_bit1_end}   &&   $other->{_bit1_st} < $other->{_bit1_end} ){
			$bit1_ovlp = overlap( $self->{_bit1_st}, $self->{_bit1_end}, $other->{_bit1_st}, $other->{_bit1_end} );
		}	
		elsif ( $self->{_bit1_st} > $self->{_bit1_end}   &&   $other->{_bit1_st} > $other->{_bit1_end} ){
			$bit1_ovlp = overlap( $self->{_bit1_end}, $self->{_bit1_st}, $other->{_bit1_end}, $other->{_bit1_st} );
		}
		else{
			return 0;
		}
		if ( $self->{_bit2_st} < $self->{_bit2_end}   &&   $other->{_bit2_st} < $other->{_bit2_end} ){
			$bit2_ovlp = overlap( $self->{_bit2_st}, $self->{_bit2_end}, $other->{_bit2_st}, $other->{_bit2_end} );
		}
		elsif ( $self->{_bit2_st} > $self->{_bit2_end}   &&   $other->{_bit2_st} > $other->{_bit2_end} ){
			$bit2_ovlp = overlap( $self->{_bit2_end}, $self->{_bit2_st}, $other->{_bit2_end}, $other->{_bit2_st} );
		}
		else{
			return 0;
		}
		if ( $bit1_ovlp > 0   &&   $bit2_ovlp > 0 ){
			return 1;
		}
	}
	return 0;
}

sub overlaps_with_gff_object {
# $hybrid -> overlaps_with_gff_object( $gff )
# returns 1 if two target sites overlap, 0 otherwise
# criteria for overlap: same gene, same position, and (either same strand or undefined strand in reference)
	my ($self, $other) = @_;
	if ( ($self->{_bit1_nm} eq $other->seq_ID()  &&  overlap( $self->{_bit1_st},$self->{_bit1_end},$other->start(),$other->end() )>0)
	||   ($self->{_bit2_nm} eq $other->seq_ID()  &&  overlap( $self->{_bit2_st},$self->{_bit2_end},$other->start(),$other->end() )>0)){
		if( !defined( $other->{_strand} ) || $other->{_strand} eq '.'  ||  $other->{_strand} eq 'NA' 
		||  ($other->{_strand} eq '+' && ($self->{_bit1_st} < $self->{_bit1_end}))
		||  ($other->{_strand} eq '-' && ($self->{_bit1_st} > $self->{_bit1_end}))){
			return 1;
		}
	}
	return 0;
}

sub get_bit_names{
	my ( $self ) = @_;
	return ( $self->{_bit1_nm}, $self->{_bit2_nm} );
}

sub get_sorted_bit_names{
	my ( $self ) = @_;
	return ( $self->{_sorted_bit_nm} );
}

sub check_hit_names {
	my ($self, $nm1, $nm2) = @_;
	if ( ( $self->{_bit1_nm} =~ $nm1 && $self->{_bit2_nm} =~ $nm2 )     ||      ( $self->{_bit1_nm} =~ $nm2 && $self->{_bit2_nm} =~ $nm1 ) ){
		return 1;
	}
	return 0;
}

sub merge_with {
# merges the current object with another one by extending the limits and increasing the count and averaging dG values
	my ($self, $other) = @_;
	$self->{_dG} = weighted_mean ( $self->{_dG}, $other->{_dG}, $self->{_count}, $other->{_count} );
	$self->{_count} += $other->{_count};

	if ( $self->{_bit1_st} < $self->{_bit1_end}   &&   $other->{_bit1_st} < $other->{_bit1_end} ){
		$self->{_bit1_st} = min ( $self->{_bit1_st}, $other->{_bit1_st} );
		$self->{_bit1_end} = max ( $self->{_bit1_end}, $other->{_bit1_end} );
	}
	elsif ( $self->{_bit1_st} > $self->{_bit1_end}   &&   $other->{_bit1_st} > $other->{_bit1_end} ){
		$self->{_bit1_st} = max ( $self->{_bit1_st}, $other->{_bit1_st} );
		$self->{_bit1_end} = min ( $self->{_bit1_end}, $other->{_bit1_end} );	
	}
	else{
		die "internal error: attempted merging features with opposite orientations";
	}
	if ( $self->{_bit2_st} < $self->{_bit2_end}   &&   $other->{_bit2_st} < $other->{_bit2_end} ){
		$self->{_bit2_st} = min ( $self->{_bit2_st}, $other->{_bit2_st} );
		$self->{_bit2_end} = max ( $self->{_bit2_end}, $other->{_bit2_end} );
	}
	elsif ( $self->{_bit2_st} > $self->{_bit2_end}   &&   $other->{_bit2_st} > $other->{_bit2_end} ){
		$self->{_bit2_st} = max ( $self->{_bit2_st}, $other->{_bit2_st} );
		$self->{_bit2_end} = min ( $self->{_bit2_end}, $other->{_bit2_end} );	
	}
	else{
		die "internal error: attempted merging features with opposite orientations";
	}

	$self->{_bit1_st_in_rd} = ".";
	$self->{_bit1_end_in_rd} = ".";
	$self->{_bit1_eval} = ".";
	$self->{_bit2_st_in_rd} = ".";
	$self->{_bit2_end_in_rd} = ".";
	$self->{_bit2_eval} = ".";

	$self->{_found_overlap}++;

	if ( $self->{_experiment} ne $other->{_experiment} ){
		$self->{_experiment} = $self->{_experiment} . "_" . $other->{_experiment};
	}

	@{ $self->{_seq_ID_list} } = (@{ $self->{_seq_ID_list} }, @{ $other->{_seq_ID_list} });

	return 1;
}

sub reverse_bit_order {
	my ($self) = @_;

	my $tmp_nm = $self->{_bit1_nm};
	$self->{_bit1_nm} = $self->{_bit2_nm};
	$self->{_bit2_nm} = $tmp_nm;

	my $tmp_st = $self->{_bit1_st};
	my $tmp_end = $self->{_bit1_end};
	$self->{_bit1_st} = $self->{_bit2_st};
	$self->{_bit1_end} = $self->{_bit2_end};
	$self->{_bit2_st} = $tmp_st;
	$self->{_bit2_end} = $tmp_end;
	
	$self->{_bit1_st_in_rd} = ".";
	$self->{_bit1_end_in_rd} = ".";
	$self->{_bit1_eval} = ".";
	$self->{_bit2_st_in_rd} = ".";
	$self->{_bit2_end_in_rd} = ".";
	$self->{_bit2_eval} = ".";

	return 1;
}

sub convert_to_plus_strand {
# if both bit1 and bit2 are in antisense orientation, then they are both converted to sense orientation
	my ($self) = @_;
	if ( ($self->{_bit1_st} > $self->{_bit1_end})  &&  ($self->{_bit2_st} > $self->{_bit2_end}) ){
		my $tmp = $self->{_bit1_nm};
		$self->{_bit1_nm} = $self->{_bit2_nm};
		$self->{_bit2_nm} = $tmp;
		my $tmp1 = $self->{_bit1_st};
		my $tmp2 = $self->{_bit1_end};
		$self->{_bit1_st} = $self->{_bit2_end};
		$self->{_bit1_end} = $self->{_bit2_st};
		$self->{_bit2_st} = $tmp2;
		$self->{_bit2_end} = $tmp1;

	        $self->{_bit1_st_in_rd} = ".";
	        $self->{_bit1_end_in_rd} = ".";
	        $self->{_bit1_eval} = ".";
	        $self->{_bit2_st_in_rd} = ".";
	        $self->{_bit2_end_in_rd} = ".";
	        $self->{_bit2_eval} = ".";
	}

	return 1;
}

sub print_coord {
	my ( $self ) = @_;
	my $line = join( "\t", $self->{_bit1_nm}, $self->{_bit1_st}, $self->{_bit1_end}, $self->{_bit2_nm}, $self->{_bit2_st}, $self->{_bit2_end});
	print $line; 
}

sub print_hyb {
	my ( $self ) = @_;
	my $dG = ".";
	if ( defined( $self->{_dG} ) ){
		$dG = sprintf ( "%.2f", $self->{_dG} );
	}
	my $line = join( "\t", $self->{_seq_ID}, $self->{_seq}, $dG, $self->{_bit1_nm}, $self->{_bit1_st_in_rd}, $self->{_bit1_end_in_rd}, $self->{_bit1_st}, $self->{_bit1_end}, $self->{_bit1_eval}, $self->{_bit2_nm}, $self->{_bit2_st_in_rd}, $self->{_bit2_end_in_rd}, $self->{_bit2_st}, $self->{_bit2_end}, $self->{_bit2_eval}, $self->{_count});
	print $line;
}

sub print_hyb_15_columns {
	my ( $self ) = @_;
	my $dG = ".";
	if ( defined( $self->{_dG} ) ){
		$dG = sprintf ( "%.2f", $self->{_dG} );
	}
	my $line = join( "\t", $self->{_seq_ID}, $self->{_seq}, $dG, $self->{_bit1_nm}, $self->{_bit1_st_in_rd}, $self->{_bit1_end_in_rd}, $self->{_bit1_st}, $self->{_bit1_end}, $self->{_bit1_eval}, $self->{_bit2_nm}, $self->{_bit2_st_in_rd}, $self->{_bit2_end_in_rd}, $self->{_bit2_st}, $self->{_bit2_end}, $self->{_bit2_eval} );
	print $line;
}

sub found_overlap {
	my ( $self, $fnd ) = @_;
	$self->{_found_overlap} = $fnd if defined( $fnd );
	return $self->{_found_overlap};
}

sub twoway_overlap {
	my ( $self, $ovlp ) = @_;
	$self->{_twoway_overlap} = $ovlp if defined( $ovlp );
	return $self->{_twoway_overlap};
}

sub seq_ID {
        my ( $self, $seq_ID ) = @_;
        $self->{_seq_ID} = $seq_ID if defined($seq_ID);
        return $self->{_seq_ID};
}

sub seq_ID_list {
	my ( $self ) = @_;	
	return @{ $self->{_seq_ID_list} };
}

sub count {
	my ( $self, $count ) = @_;
	$self->{_count} = $count if defined($count);
	return $self->{_count};
}

sub experiment {
	my ( $self, $exp ) = @_;
	$self->{_experiment} = $exp if defined($exp);
	return $self->{_experiment};
}

sub add_note {
	my ( $self, $note ) = @_;
	$self->{_note} .= $note if defined($note);
	return 1;
}

sub sum_e_values {
	my ( $self ) = @_;
	return $self->{_bit1_eval} + $self->{_bit2_eval};
}

sub match_bit_name {
	my ( $self, $nm ) = @_;
	if ( $self->{_bit1_nm} =~ $nm ){
		return $self->{_bit1_nm};
	}
        elsif ( $self->{_bit2_nm} =~ $nm ){
                return $self->{_bit2_nm};
        }
	else{
		return undef;
	}
}

sub line {
	my ( $self, $line ) = @_;
	$self->{_line} = $line if defined($line);
	return $self->{_line};
}

sub min {
	my ($a, $b) = @_;
	$a <= $b ? $a : $b;
}

sub max {
	my ($a, $b) = @_;
	$a >= $b ? $a : $b;
}

sub overlap {
	my ($a, $b, $c, $d) = @_;
	1 + min( $b, $d ) - max ($a, $c );
}

sub weighted_mean{
	my ($a, $b, $a_weight, $b_weight) = @_;
	if ( defined($a_weight)  &&  defined($b_weight) ){
		return ( ( $a * $a_weight + $b * $b_weight )   /   ( $a_weight + $b_weight) );
	}
	else{
		return ( $a + $b ) / 2;
	}
}

return 1;

