#!/usr/bin/awk -f

/^@/{
	sub(/^@/,">@",$1)
	print; getline; print; getline; getline
}

