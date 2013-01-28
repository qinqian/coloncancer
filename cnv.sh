#!/bin/bash

## preprocessing refgene to bed with symbol
#sed 1d ./processed/hg19.refgene | cut -f 3,5,6,13 > ./processed/hg19.bed

## gather data and classify simply
gather() {
    echo $0
    find $1 -name "*snp*" -type f -print0 | xargs -0 -i cp {} $2
    cd $2
    mkdir seg nocnvseg
    mv *nocnv* nocnvseg
    mv *seg.txt seg
}

## intersect to get gene to copies
seg_to_bed() {
    echo $1
    # sed 1d $1 | cut -f 2,3,4,6 | \
    # awk -F '\t' '\
    # {
    #     printf("chr%s\n", $0)
    # }' > ${1}.bed
    # intersectBed -f $2 -wa -wb  -a ${1}.bed -b ./processed/hg19.bed > ${1}genes.dat
    cut -f 4,8 ${1}genes.dat > ${1}genescopy.dat
}

## using python to get batch info
batch_average(){
    $input=$1
    $output=$2
    $ref=$3
    python cnv.py $input $ref $output
    paste $output final_cnv.txt
}

main() {
  #gather ../COAD ../snpcnv
  seg=`ls ../snpcnv/seg/*.txt`
  nocnv=`ls ../snpcnv/nocnvseg/*.txt`
  for s in $seg
  do
      seg_to_bed $s 1e-9
  done

  # batch_average 
  # time python cnv.py ../snpcnv/seg/ hg19 ../snpcnv/output/
  paste -d\\t *.txt > snp_cnv_1bp.txt
}
    
main
