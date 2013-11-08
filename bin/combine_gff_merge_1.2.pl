#!/usr/bin/perl

use warnings;
use strict;
use lib '/homes/gkudla/bin';
use Hybrid;
my @all_hits;
my $EXP = "default";
my $TARGET = ".";
my $GUIDE = ".";
my $VERBOSE = 0;
my $RESET_COUNT = 0;

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;

my $sys_sort_input = "sort -k7,7r -k1,1 -k4,4n $ARGV[0] > temp.gff\n";
system ( $sys_sort_input );
open ( IN_SORTED, ("<temp.gff") ) or die "cannot create or open file temp.gff";

LINE:	while(<IN_SORTED>){

	my $curr_hit = new Hybrid();
	$curr_hit -> initialize_gff( $_, $TARGET, $GUIDE ) or next LINE;
	if ( $curr_hit -> strand()    !~    /^[+-]$/ ){ my $strand = $curr_hit -> strand(); die "forbidden strand value $strand"; }
	if ( $RESET_COUNT ) { $curr_hit -> count( 1 ) };
	if ( $VERBOSE ) { $curr_hit -> print_simple(); print("\n"); }

##############################################################################################
# below: compare current hit with list of all previously found hits
# as of version 1.2: compare current hit with last hit only. this assumes that the input file is sorted, but is TONS faster!!!!
# the list of previously found hits is amended either by putting the current hit into it
# or by merging the current hit with one of the previously found hits
##############################################################################################

	if ( $#all_hits >= 0 ){
		my $old_hit = $all_hits[ $#all_hits ];

##############################################################################################
# the assignment above does not create a copy of the object
# Instead, $old_hit becomes a pointer to the object that is already stored in the array.
##############################################################################################		

		if ( $old_hit -> gff_overlaps( $curr_hit ) ){
			if ( $VERBOSE ) { print "OVERLAP FUNCTION:\t"; $curr_hit -> print_simple(); print "\tOVERLAPS\t"; $old_hit -> print_simple(); print "\n"; }
			$curr_hit -> found_overlap( 1 );
			$old_hit -> merge_with( $curr_hit );
			if ( $VERBOSE ) { print "MERGE FUNCTION:\t"; $old_hit -> print_simple(); print "\tCREATED\n\n"; }
		}
	}
	push ( @all_hits, $curr_hit ) unless $curr_hit -> found_overlap();
}
	
##############################################################################################
# print all the hits found
##############################################################################################

foreach my $hit( @all_hits ){
	$hit -> print_simple(); print("\n");
#	$hit -> print(); print("\n");
}

my $sys_cleanup = "rm temp.gff\n";
system ( $sys_cleanup );

