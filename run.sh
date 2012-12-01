#!/bin/bash
# main function of TCGA data process




main()
{
 Rscript TCGAcolon.R ../data/ colon_cancer_TCGA_agilent_expression.xls colon_cancer_mutation_all.maf
 #Rscript TCGAcolon.R ../data/ colon_cancer_TCGA_agilent_expression.xls colon_cancer_somatic_mutation.maf
}

main
