import os
import sys
import string

#path='/mnt/Storage/home/duanxk/TCGA/GBM/'
path='../COAD/'
dirs=os.listdir(path)
disease=path.split('/')[-2]
RNAseq={}
for d in dirs:
    if d.startswith(disease+'_Batch'):
        dirs2=os.listdir(path+d)
        if 'RNASeq' in dirs2:
            batchid=d.split('_')[1][5:]
            dirs3=os.listdir(path+d+'/RNASeq')
            for d3 in dirs3:
                if d3.startswith('UNC'):
                    #os.chdir(path+d+'/Expression-Genes/'+d3+'/Level_3')
                    dirs4=os.listdir(path+d+'/RNASeq/'+d3+'/Level_3')
                    for d4 in dirs4:
                        if d4.endswith('expression_gene.txt'):
                            patient=d4.split('__')[2][:15]
                            values=[]
                            inf=open(path+d+'/RNASeq/'+d3+'/Level_3/'+d4,'rU')
                            inf.readline()
                            for line in inf:
                                li=line.strip().split('\t')
                                if len(li)>2:
                                    #if li[2]=='null':
                                    #    print batchid
                                    #   print patient
                                    #   sys.exit(1)
                                    values.append(li[2])
                                else:
                                    pass
                            inf.close()
                            if (patient,batchid) in RNAseq.keys():
                                RNAseq[(patient,batchid)].append(values)
                            else:
                                RNAseq[(patient,batchid)]=[values]
                else:
                    pass
        else:
            pass

    else:
        pass

#################Caculate the mean of the replicates!!##############################

RNAseq_final={}

for Agi in RNAseq.keys():
    if len(RNAseq[Agi]) >1:
        value_final=[]
        for i in range(0,len(RNAseq[Agi][0])):
            sum=0
            remove=0
            for j in range(0,len(RNAseq[Agi])):
                if RNAseq[Agi][j][i]=='null':
                    remove=remove+1
    
                else:
                    sum=sum+float(RNAseq[Agi][j][i])
            value_final.append(sum/float(len(RNAseq[Agi])-remove))
        #print value_final
        RNAseq_final[Agi]=value_final
    else:
        RNAseq_final[Agi]=RNAseq[Agi][0]
RNAseq=0

######################The genes name list&&!!!Modify the sample path!!!!###########################    
RNAseq_sample='../COAD/COAD_Batch76_level3/RNASeq/UNC__IlluminaGA_RNASeq/Level_3/unc.edu__IlluminaGA_RNASeq__TCGA-A6-2671-01A-01R-1410-07__expression_gene.txt'
inf2=open(RNAseq_sample,'rU')
RNAseq_genes=[]
inf2.readline()
for line2 in inf2:
    li2=line2.split('\t')
    if len(li2)>2:
        RNAseq_genes.append(li2[1])
    else:
        pass
inf2.close()

########################Write the combat input files##########################
RNAseq_exp=open(path+'RNAseq_all_Count_expression.txt','w')
RNAseq_info=open(path+'RNAseq_sample_Count_batch_info.txt','w')
RNAseq_info.write('Array name\tSample name\tBatch\n')
Agi_keys=RNAseq_final.keys()
i=0
for key in Agi_keys:
    if i==len(Agi_keys)-1:
        RNAseq_exp.write(key[0]+'_'+key[1]+'\n')
    else:
        RNAseq_exp.write(key[0]+'_'+key[1]+'\t')
    i=i+1
    RNAseq_info.write(key[0]+'_'+key[1]+'\t'+key[0]+'_'+key[1]+'\t'+key[1]+'\n')
for i in range(0,len(RNAseq_genes)):
    RNAseq_exp.write(RNAseq_genes[i])
    
    for j in range(0,len(Agi_keys)):
        #print Agi_final
        RNAseq_exp.write('\t'+str(RNAseq_final[Agi_keys[j]][i]))
    RNAseq_exp.write('\n')

RNAseq_exp.close()
RNAseq_info.close()
