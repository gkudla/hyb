#class Gff
package Gff;
use strict;

#constructor
sub new {
	my ($class) = @_;
	my $self = {
		_seq_ID		=> undef,
		_source		=> undef,
		_feature	=> undef,
		_start		=> undef,
		_end		=> undef,
		_count		=> undef,
		_strand		=> undef,
		_frame		=> undef,
		_note		=> undef,
		_found_overlap	=> undef,
		_line		=> undef
	};
	bless $self, $class;
	return $self;
}

sub initialize_gff {
	my ( $self, $line ) = @_;
	chomp $line;
	my ( $seq_ID, $source, $feature, $start, $end, $count, $strand, $frame, $note ) = split("\t", $line, -1);
	$self->{_seq_ID} = $seq_ID if defined( $seq_ID );
	$self->{_source} = $source if defined( $source );
	$self->{_feature} = $feature if defined( $feature );
	$self->{_start} = $start if defined( $start );
	$self->{_end} = $end if defined( $end );
	$self->{_strand} = $strand if defined( $strand );
	$self->{_note} = $note if defined( $note );
	if ( defined( $count ) && $count =~ /^[0-9]+$/ ) {
		$self->{_count} = $count;
	}
	else {
		$self->{_count} = 0;
	}
	$self->{_found_overlap} = 0;

#	line not used as it uses too much memory
#	$self->{_line} = $line if defined $line; 
	return 1;
}

sub found_overlap {
	my ( $self, $fnd ) = @_;
	$self->{_found_overlap} = $fnd if defined( $fnd );
	return $self->{_found_overlap};
}

# returns 1 if two items overlap, 0 otherwise
sub overlaps {
	my ($self, $other) = @_;
	my $ovlp = overlap( $self->{_start}, $self->{_end}, $other->{_start}, $other->{_end} );
	my $strands_matched = match_strand( $self, $other );
	if ( $ovlp > 0   &&   $strands_matched > 0    &&   $self->{_seq_ID} eq $other->{_seq_ID} ){
		return 1;
	}
	else{
		return 0;
	}
}

# returns 1 if two items overlap but are on opposite strands, 0 otherwise
sub overlaps_antisense {
	my ($self, $other) = @_;
	my $ovlp = overlap( $self->{_start}, $self->{_end}, $other->{_start}, $other->{_end} );
	my $strands_matched = match_strand( $self, $other );
	if ( $ovlp > 0   &&   $strands_matched == 0    &&   $self->{_seq_ID} eq $other->{_seq_ID} ){
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
	$self->{_count} += $other->{_count};
	$self->{_start} = min ( $self->{_start}, $other->{_start} );
	$self->{_end} = max ( $self->{_end}, $other->{_end} );
	return 1;
}

sub print {
	my ( $self ) = @_;
	printf( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t.\t%s", $self->{_seq_ID}, $self->{_source}, $self->{_feature}, $self->{_start}, $self->{_end}, $self->{_count}, $self->{_strand}, $self->{_note} );
	return 1;
}

sub remove_redundancy {
	my @arr = split ( "_", $_[0] );
	my %seen;
	return join ( "_", grep !$seen{$_}++, sort( @arr )); 
}

sub seq_ID {
	my ( $self, $seq_ID ) = @_;
	$self->{_seq_ID} = $seq_ID if defined($seq_ID);
	return $self->{_seq_ID};
}

sub count_add_1 {
	my ( $self ) = @_;
	$self->{_count} = $self->{_count} + 1;
	return 1;
}

sub start {
	my ( $self ) = @_;
	return $self->{_start};
}

sub end {
	my ( $self ) = @_;
	return $self->{_end};
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

sub add_note {
	my ( $self, $note ) = @_;
	$self->{_note} .= $note if defined($note);
	return 1;
}

sub feature {
	my ( $self, $feature ) = @_;
	$self->{_feature} = $feature if defined($feature);
	return $self->{_feature};
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
	return undef if $self->{_seq_ID} ne $other->{_seq_ID};
        my $dist = 1 - overlap( $self->{_start}, $self->{_end}, $other->{_start}, $other->{_end} );
	return 0 if $dist <= 0;
	return -$dist if $other->{_start} < $self->{_start};
	return $dist;
}

return 1;

