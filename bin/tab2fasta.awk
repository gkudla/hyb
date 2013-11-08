#!/usr/bin/awk -f

BEGIN{
	FS = "\t"
}

!/^#/ && /^>/{print $1 "\n" $2; next}
!/^#/ && !/^>/{print ">" $1 "\n" $2}
