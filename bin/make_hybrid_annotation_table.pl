#!/usr/bin/perl

# example: make_hybrid_annotation_table.pl in.1 in.2 > out.txt
# reads any number of files in the format:
# seq_ID TAB property1=xxx;property2=yyy;...
# outputs table with seq_ID's in rows and properties in columns, for all files for all seq_ID's

$SORT_COLUMNS = 0;
$OUTPUT_ZERO_FOR_NA = 0;

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;

while (<>){

	chomp;
	my ($name,$properties) = split("\t", $_, -1);
	@curr_hyb_properties = split(";", $properties, -1);
	
	my $a = undef;
	if ( !($all_hyb_names{$name}) ){
		$all_hyb_names{$name} = 1;
		push @all_hyb_names, $name;
	}
	while (($a = shift @curr_hyb_properties) && $a =~ /([A-Za-z_0-9-+]+)=([^;]+)/){
#	while (($a = shift @curr_hyb_properties) && $a =~ /([A-Za-z_0-9]+)=([^;]+)/){
	        if ( !($all_hyb_properties{$1}) ){
	                $all_hyb_properties{$1} = 1;
			push @all_hyb_properties, $1;
		}
		$data{$name}{$1} = $2;
	}
}

if( $SORT_COLUMNS ){
	@all_hyb_properties = sort @all_hyb_properties;
}

##################################
# print header line
##################################

print "seq_ID";
foreach my $property (@all_hyb_properties){
	print "\t$property";
}
print "\n";

##################################
# print data
##################################

foreach my $name (@all_hyb_names){
	print "$name";
	foreach my $property (@all_hyb_properties){
		my $curr_txt = $data{$name}{$property};
		if( defined $data{$name}{$property} ){
			$curr_txt = $data{$name}{$property}
		}
		else{	
			$curr_txt = "0" if ($OUTPUT_ZERO_FOR_NA==1);
			$curr_txt = "NA" if ($OUTPUT_ZERO_FOR_NA==0);
		}
		print "\t$curr_txt";
#		print "\t$data{$name}{$property}";
	}
	print "\n";
}	


