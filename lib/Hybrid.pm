#class Hybrid
package Hybrid;
use strict;

#constructor
sub new {
	my ($class) = @_;
	my $self = {
		_experiment	=> undef,
		_target		=> undef,
		_guide		=> undef,
		_start		=> undef,
		_end		=> undef,
		_count		=> undef,
		_strand		=> undef,
		_dG		=> undef,
		_found_overlap	=> undef,
		_note		=> undef,
		_line		=> undef
	};
	bless $self, $class;
	return $self;
}

#accessor method for Hybrid elements
sub initialize_hyb {
	my ( $self, $line, $exp, $target, $guide ) = @_;
	my ( $seq_ID, $seq, $dG, $bit1_nm, $bit1_st_in_rd, $bit1_end_in_rd, $bit1_st, $bit1_end, $bit1_eval, $bit2_nm, $bit2_st_in_rd, $bit2_end_in_rd, $bit2_st, $bit2_end, $bit2_eval, $count) = split("\t", $line, -1);
	if ( $bit1_nm =~ $target && $bit2_nm =~ $guide ){
		$self->{_target} = $bit1_nm if defined( $bit1_nm );
		$self->{_guide} = $bit2_nm if defined( $bit2_nm );
		$self->{_start} = $bit1_st if defined( $bit1_st );
		$self->{_end} = $bit1_end if defined( $bit1_end );
	}
	elsif ( $bit1_nm =~ $guide && $bit2_nm =~ $target ){
		$self->{_target} = $bit2_nm if defined( $bit2_nm );
		$self->{_guide} = $bit1_nm if defined( $bit1_nm );
		$self->{_start} = $bit2_st if defined( $bit2_st );
		$self->{_end} = $bit2_end if defined( $bit2_end );
	}
	else {
		return 0;
	}
	$self->{_line} = $line if defined( $line );
	$self->{_experiment} = $exp if defined( $exp );
	$self->{_dG} = $dG if defined( $dG );
	if ( defined( $count ) && $count =~ /^[0-9]+$/ ) {
		$self->{_count} = $count;
	}
	else {
		$self->{_count} = 1;
	}
	$self->{_found_overlap} = 0;
	return 1;
}

sub initialize_gff {
	my ( $self, $line, $target_to_match, $guide_to_match ) = @_;
	chomp $line;
	my ( $target, $experiment, $guide, $start, $end, $count, $strand, $frame, $notes ) = split("\t", $line, -1);
	if ( $target =~ $target_to_match && $guide =~ $guide_to_match ){
		$self->{_line} = $line if defined( $line );
		$self->{_experiment} = $experiment if defined( $experiment );
		$self->{_target} = $target if defined( $target );
		$self->{_guide} = $guide if defined( $guide );
		$self->{_start} = $start if defined( $start );
		$self->{_end} = $end if defined( $end );
		$self->{_strand} = $strand if defined( $strand );
		$self->{_note} = $notes if defined( $notes );
		if ( defined( $count ) && $count =~ /[0-9]+/ ) {
			$self->{_count} = $count;
		}
		else {
			$self->{_count} = 0;
		}
		if ( defined( $notes ) && $notes =~ /dG=([0-9.-]+);/ ){
			$self->{_dG} = $1;
		}
		$self->{_found_overlap} = 0;
	}
	else {
		return 0;
	}
	return 1;
}

sub found_overlap {
	my ( $self, $fnd ) = @_;
	$self->{_found_overlap} = $fnd if defined( $fnd );
	return $self->{_found_overlap};
}

# returns 1 if two target sites overlap, 0 otherwise
sub overlaps {
	my ($self, $other) = @_;
	my $ovlp = overlap( $self->{_start}, $self->{_end}, $other->{_start}, $other->{_end} );
	my $strands_matched = match_strand( $self, $other );
	if ( $ovlp > 0   &&   $strands_matched > 0    &&   $self->{_target} eq $other->{_target}   &&   $self->{_guide} eq $other->{_guide} ){
		return 1;
	}
	else{
		return 0;
	}
}

# changed overlap criteria: match in guide (3rd field of gff file) not required. returns 1 if two target sites overlap, 0 otherwise
sub gff_overlaps {
	my ($self, $other) = @_;
	my $ovlp = overlap( $self->{_start}, $self->{_end}, $other->{_start}, $other->{_end} );
	my $strands_matched = match_strand( $self, $other );
	if ( $ovlp > 0   &&   $strands_matched > 0    &&   $self->{_target} eq $other->{_target} ){
#	if ( $ovlp > 0   &&   $self->{_target} eq $other->{_target}   &&   ( !defined($self->{_strand}) || !defined($other->{_strand}) || $self->{_strand} eq $other->{_strand} ) ){
		return 1;
	}
	else{
		return 0;
	}
}

sub match_strand {
	my ($self, $other) = @_;
	if ( !defined( $self->{_strand} ) || !defined( $other->{_strand} ) ){
		return 1;
	}	
	if ( $self->{_strand} eq '.'  || $other->{_strand} eq '.'  ||  $self->{_strand} eq 'NA'  || $other->{_strand} eq 'NA' ){
		return 1;
	}	
	if ( $self->{_strand} eq $other->{_strand} ){
		return 1;
	}	
	return 0;	
}

# merges the current object with another one by extending the limits and increasing the count and dG values
sub merge_with {
	my ($self, $other) = @_;
	$self->{_dG} = weighted_mean ( $self->{_dG}, $other->{_dG}, $self->{_count}, $other->{_count} );
	$self->{_count} += $other->{_count};
	$self->{_start} = min ( $self->{_start}, $other->{_start} );
	$self->{_end} = max ( $self->{_end}, $other->{_end} );
	if ( $self->{_experiment} ne $other->{_experiment} ){
		$self->{_experiment} = $self->{_experiment} . "_" . $other->{_experiment};
	}
	return 1;
}

sub print {
	my ( $self ) = @_;
	my $note;
	my $dG;
	if ( defined( $self->{_note}) && defined( $self->{_dG} ) ){
		$self->{_note} =~ s/([^;])$/$1;/;
		$self->{_note} =~ s/dG=[0-9.-]+;//;
		$dG = sprintf ( "dG=%.2f;", $self->{_dG} );
		$note = $dG . $self->{_note} ;
	}
	elsif ( defined( $self->{_note}) && !defined( $self->{_dG} ) ){
		$self->{_note} =~ s/([^;])$/$1;/;
		$note = $self->{_note} ;
	}
	elsif ( !defined( $self->{_note}) && defined( $self->{_dG} ) ){
		$dG = sprintf ( "dG=%.2f;", $self->{_dG} );
		$note = $dG ;
	}
	else{
		$note = ".";
	}
	printf( "%s\t%s\t%s\t%s\t%s\t%s\t.\t.\t%s", $self->{_target}, $self->{_experiment}, $self->{_guide}, $self->{_start}, $self->{_end}, $self->{_count}, $note );
	return 1;
}

sub print_simple {
	my ( $self ) = @_;
	my $experiment = &remove_redundancy( $self->{_experiment} );
	printf( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t.", $self->{_target}, $experiment, $self->{_guide}, $self->{_start}, $self->{_end}, $self->{_count}, $self->{_strand} );
	return 1;
}

sub print_simple_with_note {
	my ( $self ) = @_;
	my $experiment = &remove_redundancy( $self->{_experiment} );
	printf( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t%s", $self->{_target}, $experiment, $self->{_guide}, $self->{_start}, $self->{_end}, $self->{_count}, $self->{_strand}, $self->{_note} );
	return 1;
}

sub remove_redundancy {
	my @arr = split ( "_", $_[0] );
	my %seen;
	return join ( "_", grep !$seen{$_}++, sort( @arr )); 
}

sub target {
	my ( $self, $target ) = @_;
	$self->{_target} = $target if defined($target);
	return $self->{_target};
}

sub count_add_1 {
	my ( $self ) = @_;
	$self->{_count} = $self->{_count} + 1;
	return 1;
}

sub count {
	my ( $self, $count ) = @_;
	$self->{_count} = $count if defined($count);
	return $self->{_count};
}

sub strand {
	my ( $self, $strand ) = @_;
	$self->{_strand} = $strand if defined($strand);
	return $self->{_strand};
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

sub gff_distance{
# calculates distance between two features on a chromosome, disregarding strand
# returns undef if the features are on different chromosomes
# returns 0 if the features overlap
# returns a negative number, the distance between features, if the $other feature is upstream
# returns a positive number if the $other feature is downstream
# if the two features are next to each other, but don't overlap, the distance is one 
        my ($self, $other) = @_;
	return undef if $self->{_target} ne $other->{_target};
        my $dist = 1 - overlap( $self->{_start}, $self->{_end}, $other->{_start}, $other->{_end} );
	return 0 if $dist <= 0;
	return -$dist if $other->{_start} < $self->{_start};
	return $dist;
}

return 1;

