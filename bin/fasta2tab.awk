#!/usr/bin/awk -f

BEGIN{
	first_line=1
}

/^#/{
	next	
}

$0==""{
	next
}

/^>/ && first_line {
	sub(/>/,"",$0)
	printf "%s\t", $0
	first_line=0
	next
}

/^>/ {
	sub(/>/,"",$0)
	printf "\n%s\t", $0
	next
}

{
	printf "%s", $0
}

END{
	printf "\n"
}
