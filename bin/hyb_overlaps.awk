#!/usr/bin/awk -f

# This script reads a file in .hyb format, and calculates the overlap (L) between the arms of homotypic chimeras, as described in Gabryelska et al (2022)
# The scripts outputs a list of homotypic chimeras with the value of L printed in the final column (column 16)

# Usage: hyb_overlaps.awk input.hyb > output.hyb.overlaps


BEGIN{ OFS="\t" }

$4==$10{
        ovlp=overlap( $7, $8, $13, $14)
        print $0 "\t" ovlp
}

function min(a, b){
        if(a<b){return a}
        return b
}

function max(a, b){
        if(a>b){return a}
        return b
}

function overlap(a, b, c, d){
        return 1 + min( b, d ) - max( a, c )
}
