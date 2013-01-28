#!/bin/bash
# not work for array in row pattern
# by cut and merge by step
# cat test1 > outputfile;for file in `ls test* `;do cut -d" " -f2 $file > temp; paste -d" " outputfile temp > outputfile;done
## the same
#temp=$(cat aaa)
#echo $temp
#for i in test*; 
#do
#    echo aaa | join -j 1 - $i
#done;

#join -j 1 test1 test2  | join -j 1 test3 - | join -j 1 test4 -
## by cut together and paste once

## ways to merge on common fields


i=0
sed 1d snp_cnvTCGA-AA-3556-01.txt | cut -f 1 | sed 1i\
"genes" > delim
for file in snp*TCGA*
do
	i=$(($i+1))
	cut -f 2 $file > ${file}__${i}.temp
done
paste -d\\t delim snp_cnv*__*.temp > output
