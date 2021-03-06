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
            printf("%s\t%s\t%s\t%s\tExon_%d\n",$1,$13,$3, int(exonEnds[i])-int(exonStarts[i]), ($3=="+"?i:exonCount-i))
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

main()
{
    #Rscript -e ""
    # Rscript TCGAcolon.R ../data/
    #colon_cancer_TCGA_agilent_expression.xls
    #colon_cancer_mutation_all.maf

    # getExon ../ref/hg19_allfields.txt ${temp}/hg19_exon.ref


    ## for agilent expression data remove null and NULL


    ## for RNASeq, remove |ID ?? and ? gene
    sed -e '/\?/d' -e 's/|[0-9]*//g' ../COAD/RNAseq_all_expression.txt > ../data/RNAseq_all_expression_Trim.txt

    ## remove replicates gene names
    cut -f 1 ../data/RNAseq_all_expression_Trim.txt | sort -d | uniq > RNAsequniq
    cut -f 1 ../data/RNAseq_all_expression_Trim.txt > RNAseqName
    diff RNAseqName RNAsequniq ## gene SLC35E2 has replicates, named
                               ## as _1 and _2,  focus on this gene if possible

    #sed 's/\([[:alpha:]]\)\([^ \n]*\)/\2\1ay/g' testsed
    #/black/!s/cow/horse/ testsed

    # Rscript TCGAcolon.R ../data/ colon_cancer_TCGA_agilent_expression.xls colon_cancer_somatic_mutation.maf
    #Rscript TCGAcolon.R ../data/ colon_cancer_TCGA_agilent_expression.xls colon_cancer_mutation_all.maf
}
main

