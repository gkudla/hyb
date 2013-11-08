# input hybrid.hit_table.txt

echo $1

echo "all	"	`awk '{a+=$3}END{print a}' $1`

echo "rRNA-snR*	"	`awk '($1~/RDN37/&&$2~/snR/)||($2~/RDN37/&&$1~/snR/){a+=$3}END{print a}' $1` 

echo "other_A-B	"	`awk '$1!=$2&&($1!~/RDN37/||$2!~/snR/){a+=$3}END{print a}' $1`

echo "rRNA-rRNA	"	`awk '($1~/RDN37/&&$2~/RDN37/&&$1==$2){a+=$3}END{print a}' $1`

echo "snR*-snR*	"	`awk '($1~/snR/&&$2~/snR/&&$1==$2){a+=$3}END{print a}' $1`

echo "other_A-A	"	`awk '$1==$2&&$1!~/RDN37/&&$2!~/snR/{a+=$3}END{print a}' $1`

echo

awk '/LSR1|snR6[^0-9]|snR14/' $1






 
