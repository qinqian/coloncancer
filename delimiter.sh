IN="bla@some.com;john@home.com"
arr=$(echo $IN | tr ";" "\n")

## change default delimiter
OIFS=$IFS
IFS=';'
arra2=$IN

## recover
IFS=$OIFS

arr=(${str//,/})