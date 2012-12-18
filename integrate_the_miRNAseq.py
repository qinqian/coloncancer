import os
import sys
import string

#path='/mnt/Storage/home/duanxk/TCGA/GBM/'
path='/mnt/Storage/home/duanxk/TCGA/BRCA/'
dirs=os.listdir(path)
disease=path.split('/')[-2]
miRNA={}
for d in dirs:
    if d.startswith(disease+'_Batch'):
        dirs2=os.listdir(path+d)
        if 'RNASeq' in dirs2:
            batchid=d.split('_')[1][5:]
            dirs3=os.listdir(path+d+'/RNASeq')
            for d3 in dirs3:
                if d3.startswith('BCGSC'):
                    #os.chdir(path+d+'/Expression-Genes/'+d3+'/Level_3')
                    dirs4=os.listdir(path+d+'/miRNASeq/'+d3+'/Level_3')
                    for d4 in dirs4:
                        if d4.endswith('mirna_quantification.txt'):
                            patient=d4.split('__')[2][:15]
                            values=[]
                            inf=open(path+d+'/miRNASeq/'+d3+'/Level_3/'+d4,'rU')
                            inf.readline()
                            for line in inf:
                                li=line.strip().split('\t')
                                if len(li)>2:
                                    #if li[2]=='null':
                                    #    print batchid
                                    #   print patient
                                    #   sys.exit(1)
                                    values.append(li[3])
                                else:
                                    pass
                            inf.close()
                            if (patient,batchid) in miRNA.keys():
                                miRNA[(patient,batchid)].append(values)
                            else:
                                miRNA[(patient,batchid)]=[values]
                else:
                    pass
        else:
            pass

    else:
        pass

#################Caculate the mean of the replicates!!##############################

miRNA_final={}

for Agi in miRNA.keys():
    if len(miRNA[Agi]) >1:
        value_final=[]
        for i in range(0,len(miRNA[Agi][0])):
            sum=0
            remove=0
            for j in range(0,len(miRNA[Agi])):
                if miRNA[Agi][j][i]=='null':
                    remove=remove+1
    
                else:
                    sum=sum+float(miRNA[Agi][j][i])
                
            value_final.append(sum/float(len(miRNA[Agi])-remove))
        #print value_final
        miRNA_final[Agi]=value_final
    else:
        miRNA_final[Agi]=miRNA[Agi][0]
miRNA=0

######################The genes name list&&!!!Modify the sample path!!!!###########################    
miRNA_sample='/mnt/Storage/home/duanxk/TCGA/BRCA/BRCA_Batch109_level3/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/bcgsc.ca__IlluminaHiSeq_miRNASeq__TCGA-EW-A1J6-01A-11R-A13P-13__mirna_quantification.txt'
inf2=open(miRNA_sample,'rU')
miRNA_genes=[]
inf2.readline()
for line2 in inf2:
    li2=line2.split('\t')
    if len(li2)>2:
        miRNA_genes.append(li2[1])
    else:
        pass
inf2.close()

########################Write the combat input files##########################
miRNA_exp=open(path+'miRNA_all_expression.txt','w')
miRNA_info=open(path+'miRNA_sample_batch_info.txt','w')
miRNA_info.write('Array name\tSample name\tBatch\n')
Agi_keys=miRNA_final.keys()
i=0
for key in Agi_keys:
    if i==len(Agi_keys)-1:
        miRNA_exp.write(key[0]+'_'+key[1]+'\n')
    else:
        miRNA_exp.write(key[0]+'_'+key[1]+'\t')
    i=i+1
    miRNA_info.write(key[0]+'_'+key[1]+'\t'+key[0]+'_'+key[1]+'\t'+key[1]+'\n')
for i in range(0,len(miRNA_genes)):
    miRNA_exp.write(miRNA_genes[i])
    
    for j in range(0,len(Agi_keys)):
        #print Agi_final
        miRNA_exp.write('\t'+str(miRNA_final[Agi_keys[j]][i]))
    miRNA_exp.write('\n')

miRNA_exp.close()
miRNA_info.close()


                    
                    
                            
