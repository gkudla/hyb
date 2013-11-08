#!/usr/bin/awk
# @(#)analyse_miRNA_folding_viennad_GK.awk  2013-05-15  Last modified by A.J.Travis

# by GK
# usage: awk --posix -f /storage/home/gkudla/bin/analyse_miRNA_folding_viennad_GK.awk input.viennad > output.viennad_analysis

BEGIN{
	VERBOSE=0
	REPORT_MORE_FEATURES=0
	GET_ALL_TARGETS_NOT_JUST_mRNA=0
}

{
# processing hybrid name
	if (VERBOSE) {print}
	split ($0, a, "_")
	seq_ID = (a[1] "_" a[2])
	getline

# processing entire sequence
	entire_seq = $1
	getline

# processing microRNA sequence
#	print
	if ($0 ~ /_microRNA/){
		match ($1, /[A-Z]+/)
		pos_miRNA_start = RSTART
		miRNA_length = RLENGTH
		miRNA_seq = substr($1, RSTART, RLENGTH)
		miRNA_seed = substr ($1, 2, 6)
	}
	getline

# processing mRNA sequence
        if ($0 ~ /_mRNA/ || GET_ALL_TARGETS_NOT_JUST_mRNA){
                match ($1, /[A-Z]+/)
                pos_mRNA_start = RSTART
                mRNA_length = RLENGTH
                mRNA_seq = substr($1, RSTART, RLENGTH)
        }
	getline

# processing folding information
	if (VERBOSE) {print}
	folding_energy = $2
	gsub(/[)(]/, "", folding_energy)
	if ($1 ~ /^[().]+$/){
		miRNA_fold = substr($1, pos_miRNA_start, miRNA_length)
		seed_fold = substr($1, 2, 6)
		shifted_seed_fold = substr($1, 3, 6)
	}
	pos_first_bracket = index($1, "(")
		match($1,"[)][.]*$")
	pos_last_bracket = RSTART
	if (pos_first_bracket){
		pos_opposite_start_miRNA = pos_first_bracket+pos_last_bracket-1
		nucl_opposite_start_miRNA = substr(entire_seq,pos_opposite_start_miRNA,1)
		if (VERBOSE) {print entire_seq}
		if (VERBOSE) {print "OPPOSITE_START_miRNA:", pos_first_bracket, pos_last_bracket, pos_opposite_start_miRNA, nucl_opposite_start_miRNA, miRNA_seed, seed_fold}
		cnt_nucl_opposite_start_miRNA[nucl_opposite_start_miRNA]++
	}	
	getline

# printing miRNA sequence and fold
	if (VERBOSE) {print  miRNA_seq "\n" miRNA_fold}

# -------------------checking type of target site----------------------

	temp1 = miRNA_fold
	gsub ( /[^(]/, "", temp1 )
	num_brackets = length( temp1 )
	
 	temp2 = seed_fold
        gsub ( /[.]/, "", temp2 )
        seed_brackets = length( temp2 )

 	temp3 = shifted_seed_fold
        gsub ( /[.]/, "", temp3 )
        shifted_seed_brackets = length( temp3 )

# seed matches
	
	if (nucl_opposite_start_miRNA == "A" && miRNA_fold ~ /^.[(]{7}/ ){
                seed_type="S8"
	}
	else if (nucl_opposite_start_miRNA == "A" && miRNA_fold ~ /^.[(]{6}/){
		seed_type="S7A"
	}
	else if (miRNA_fold ~ /^.[(]{7}/ ){
		seed_type="S7"
	}
	else if (miRNA_fold ~ /^.[(]{6}/){
		seed_type="S6"
	}
	else if (miRNA_fold ~ /^..[(]{6}/){
		seed_type="S6S"
	}	
	else if (miRNA_fold ~ /^...[(]{11}|^....[(]{11}/){
		seed_type="CeS"
	}
	else{
		seed_type="none"
	}

# longest stem anywhere in microRNA
	longest_stem=0
	while ( match( miRNA_fold, /[(]+/ ) ){
		if(RLENGTH > longest_stem){
			longest_stem=RLENGTH
		}
		miRNA_fold=substr( miRNA_fold, RSTART+RLENGTH )
	}


##############################################################################################################

	if( REPORT_MORE_FEATURES == 0 ){
		print (seq_ID "\tnum_basepairs=" num_brackets ";seed_basepairs=" seed_brackets ";seed_type=" seed_type ";folding_energy=" folding_energy)
	}
	else{
		print (seq_ID "\tnum_basepairs=" num_brackets ";seed_basepairs=" seed_brackets ";shifted_seed_basepairs=" shifted_seed_brackets ";seed_type=" seed_type ";folding_energy=" folding_energy ";longest_stem=" longest_stem)
	}
	if (VERBOSE){ print ""}
}
 
END{
}
