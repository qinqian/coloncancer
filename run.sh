#!/bin/bash
# main function of TCGA data process
result="../results/"
temp="../temp/"

getExon()
{
    # get exons annotation from UCSC hg19_allfields.txt
    # $1 input annotation, $2 output file name
    # cut -f 9,10,11,13 $1 > $2
    awk -F '\t' '\
    {
	exonCount=int($9); split($10,exonStarts,"[,]"); split($11,exonEnds,"[,]");
	for(i=1;i<exonCount;i++)
	{
            exonlen[$13] += int(exonEnds[i])-int(exonStarts[i])
            #printf("%s\t%s\t%s\t%s\tExon_%d\n",$1,$13,$3, int(exonEnds[i])-int(exonStarts[i]), ($3=="+"?i:exonCount-i))
	}
    }
    END{
    for(i in exonlen)
        {
            printf("%s\t%s\n", i, exonlen[i])
        }
    }' $1 \
    >  $2
}

ExonLength()
{
    echo
}

main()
{
    #Rscript -e ""
    # Rscript TCGAcolon.R ../data/
    #colon_cancer_TCGA_agilent_expression.xls
    #colon_cancer_mutation_all.maf

    # getExon ../hg19_allfields.txt ${temp}/hg19_exon.ref

    Rscript TCGAcolon.R ../data/ colon_cancer_TCGA_agilent_expression.xls colon_cancer_somatic_mutation.maf
}
main

