#!/bin/bash

# input: .viennad file 
# USED BY: process_mim_hybrids_ds.sh

fname=$1
barename=${fname%.*}
extension=${fname##*.}



INPUT=$fname 
OUTPUT_ANALYSIS=$INPUT"_analysis"
OUTPUT_REPORT=$INPUT"_report"


awk --posix -f /homes/olahelwak/bin/analyse_miRNA_folding_viennad.awk $INPUT > $OUTPUT_ANALYSIS

echo -e "input file = "$fname '\n'  > $OUTPUT_REPORT

num_hybrids=`cat $OUTPUT_ANALYSIS | grep -c microRNA`
echo "number of analysed hybrids =" $num_hybrids >> $OUTPUT_REPORT 

num_S8_sites=`cat $OUTPUT_ANALYSIS | grep -c S8`
echo "number of S8 sites = "$num_S8_sites >> $OUTPUT_REPORT

num_S7A_sites=`cat $OUTPUT_ANALYSIS | grep -c S7A`
echo "number of S7A sites = "$num_S7A_sites >> $OUTPUT_REPORT

num_S7_sites=`cat $OUTPUT_ANALYSIS | grep -c S7`
echo "number of S7 sites = "$num_S7_sites >> $OUTPUT_REPORT

num_S6_sites=`cat $OUTPUT_ANALYSIS | grep -c S6`
echo "number of S6 sites = "$num_S6_sites >> $OUTPUT_REPORT

num_S6S_sites=`cat $OUTPUT_ANALYSIS | grep -c S6S`
echo "number of S6S sites = "$num_S6S_sites >> $OUTPUT_REPORT

num_CeS_sites=`cat $OUTPUT_ANALYSIS | grep -c CeS`
echo "number of CeS sites = "$num_CeS_sites >> $OUTPUT_REPORT

num_ClS_sites=`cat $OUTPUT_ANALYSIS | grep -c ClS`
echo "number of ClS sites = "$num_ClS_sites >> $OUTPUT_REPORT

num_5S_sites=`cat $OUTPUT_ANALYSIS | grep -c 5S`
echo "number of chimeras with minimum 5 bases paired = "$num_5S_sites >> $OUTPUT_REPORT

