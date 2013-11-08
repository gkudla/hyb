#!/usr/bin/perl

# bug fix 20131019: chimeras with too-large gap or overlap are now filtered correctly

# bug fix 20131015: chimeras now printed for last read_ID in the file

# bug found by Jon on 20130515, not solved at present: when there are three non-overlapping hits in the read, chimeras are sometimes called between the first and third bits.

# bug fix 20130412: changed variable name in print_hyb_format from @Fld to @F, this solves problem where the @Fld variable from the main function was overwritteni in print_hyb_format. 

# NOTES 20130128:
# chimeras are reported for those reads where at most 10 hits pass the MAX_OVERLAP and BLAST_THRESHOLD cutoffs - not for those reads with at most 10 hits overall
# bug: no chimeras are printed out for the last read_ID in the file

# change as of 20110321: the program now rejects all antisense hits in MODE==2 and MODE==3, independently of e-value or overlap

# .hyb format:
# SEQ_ID        SEQ     dG      BIT1_ID BIT1_st_rd      BIT1_en_rd      BIT1_st_gene    BIT1_en_gene    BIT1_eval       BIT2_ID BIT2_st_rd      BIT2_en_rd      BIT2_st_gene    BIT2_en_gene    BIT2_eval

# MODE = 0 allows antisense hits, removes hits with too much overlap ( use for nonencoded oligoA tails )
# $MODE = 1 allows antisense hits, removes hits with too much overlap or too long gaps ( use for CLASH when mapping to genomic sequence )
# $MODE = 2 allows sense hits only, removes hits with too much overlap or too long gaps ( use for CLASH when mapping to transcripts )
# $MODE = 3 allows sense hits only, removes hits with too much overlap
 
$MODE = 2;
$MAX_OVERLAP = 4;
$BLAST_THRESHOLD = 0.001;
$MAX_HITS_PER_SEQUENCE = 10;
$OUTPUT_FORMAT = "HYB";
$FS = "\t";
$, = "\t";

eval '$'.$1.$2 while ($ARGV[0] =~ /^([A-Za-z_0-9]+=)([A-Za-z\d.]+)/) && shift; # process any FOO=bar switches

print "#$0\n";
foreach $argnum (0 .. $#ARGV) { print "#$ARGV[$argnum]\n"; }
print "#MODE = $MODE\n#MAX_OVERLAP = $MAX_OVERLAP\n#BLAST_THRESHOLD = $BLAST_THRESHOLD\n";
print "#MAX_HITS_PER_SEQUENCE = $MAX_HITS_PER_SEQUENCE\n#OUTPUT_FORMAT = $OUTPUT_FORMAT\n\n";

while (<>) {
	
	@Fld = split /$FS/,$_;

# if current line has different Solexa_seq_ID from previous line, reset Solexa-seq_ID-specific variables
# if more than 1 hybrid bits were found, and at most MAX_HITS_PER_SEQUENCE of hits were found, print the lines

		if ($Fld[0] ne $name){
	      		if ($nbits > 1 && scalar @arr <= $MAX_HITS_PER_SEQUENCE ) {
				@sorted = 	map {$_->[0]}
						sort {$a->[7] <=> $b->[7] }
						map {chomp;[$_,split(/[ \t]+/)]} @arr;
				if ( $OUTPUT_FORMAT eq "BLAST" ){
					&print_blast_format( @sorted );
				}
				elsif ( $OUTPUT_FORMAT eq "HYB" ){
					&print_hyb_format( @sorted );
				}
			}	   

			@arr = ();
			$name = $Fld[0];
			$min_start=$Fld[6];
			$max_end=$Fld[7];
			$curr_overlap = 0;
			$curr_blast_threshold = 0;
			$nbits = 0;
		}

# if current line has same Solexa_seq_ID as previous line, calculate overlap and update start and end of mapped region of read
		else{
			$curr_overlap = overlap($min_start,$max_end,$Fld[6],$Fld[7]);
			if( $Fld[10] <= $BLAST_THRESHOLD ){
				$min_start = min($min_start,$Fld[6]);
				$max_end = max($max_end,$Fld[7]);
			}
		}

# debugging
#		print $_;
#		print "EVALUE", $Fld[10], $BLAST_THRESHOLD, "ORIENTATION", $Fld[8], $Fld[9], "OVERLAP", $curr_overlap;
#		print "\n";

# do not allow antisense hits if MODE==2 or MODE==3
		if( $Fld[8]>$Fld[9] && ( $MODE == 2 || $MODE ==3 )){ next; }

# remember lines with blast score at least as good as the best blast score 
		if ( $Fld[10] <= $curr_blast_threshold ){
			push @arr, $_;
		}

# skip lines with too much overlap, too low blast score, or antisense-oriented hits
		if ( $Fld[10] > $BLAST_THRESHOLD ){ next; }
		if ( $MODE == 0 && ( $curr_overlap > $MAX_OVERLAP ) ){ next; }
		if ( $MODE == 1 && ( $curr_overlap > $MAX_OVERLAP || $curr_overlap < -($MAX_OVERLAP) ) ){ next; }
		if ( $MODE == 2 && ( $curr_overlap > $MAX_OVERLAP || $curr_overlap < -($MAX_OVERLAP) ) ){ next; }
		if ( $MODE == 3 && ( $curr_overlap > $MAX_OVERLAP ) ){ next; }

# if criteria met (if not skipped), increase nbits (number of bits in the current hybrid)
# add current line to array (unless it was already added per curr_blast_threshold)
		$nbits++;
		if ( $Fld[10] != $curr_blast_threshold ){
			$curr_blast_threshold = $Fld[10];
			push @arr, $_
		}
}

if (eof && $nbits > 1 && scalar @arr <= $MAX_HITS_PER_SEQUENCE ) {
	@sorted = 	map {$_->[0]}
			sort {$a->[7] <=> $b->[7] }
			map {chomp;[$_,split(/[ \t]+/)]} @arr;
	if ( $OUTPUT_FORMAT eq "BLAST" ){
		&print_blast_format( @sorted );
	}
	elsif ( $OUTPUT_FORMAT eq "HYB" ){
		&print_hyb_format( @sorted );
	}
}
	   


sub ratio {
	my($a, $b) = @_;
	if ($b == 0) {
	return 0;
	}
	$a / $b;
}

sub min{
	my($a, $b) = @_;
	$a <= $b ? $a : $b;
}

sub max{
	my($a, $b) = @_;
	$a >= $b ? $a : $b;
}

sub overlap{
	my($a, $b, $c, $d) = @_;
	1 + min( $b, $d ) - max ($a, $c );
}

sub print_blast_format{
	my $num_lines = scalar @_ - 1;
	foreach $argnum ( 0 .. $num_lines ){
		print "@_[$argnum]\n";	
	}
	print "\n";
}

sub print_hyb_format{

# puts all the data into array a[ blast_line_index ][ blast_field_index ]
	my $num_lines = scalar @_ - 1;
	my @all_lines = @_;
	foreach $i ( 0 .. $num_lines ){
		@F = split /$FS/, $all_lines[ $i ];
		foreach $j ( 0 .. 10 ){
			$a[ $i ][ $j ] = $F[ $j ];
		}
	}

# for each potential hybrid, check overlap between bits and prints hybrids with overlap <= MAX_OVERLAP
	foreach $i ( 0 .. $num_lines ){
		foreach $j ( $i+1 .. $num_lines ){
                        $abs_curr_overlap = abs overlap( $a[$i][6],$a[$i][7],$a[$j][6],$a[$j][7]);
			if ( $abs_curr_overlap <= $MAX_OVERLAP ){	
				print "$a[$i][0]\t.\t.\t$a[$i][1]\t$a[$i][6]\t$a[$i][7]\t$a[$i][8]\t$a[$i][9]\t$a[$i][10]\t";
				print                  "$a[$j][1]\t$a[$j][6]\t$a[$j][7]\t$a[$j][8]\t$a[$j][9]\t$a[$j][10]\n";
			}
		}
	}
}
