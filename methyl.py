## methylation
import os
import sys
import string

path='../COAD/'
dirs=os.listdir(path)
disease=path.split('/')[-2]
methyl={}
for d in dirs:
    if d.startswith(disease+'_Batch'):
        dirs2=os.listdir(path+d)
        if 'DNA_Methylation' in dirs2:
            batchid=d.split('_')[1][5:]
            dirs3=os.listdir(path+d+'/DNA_Methylation')
            for d3 in dirs3:
                if d3.startswith('JHU'):
                    dirs4=os.listdir(path+d+'/DNA_Methylation/'+d3+'/Level_3')
                    for d4 in dirs4:
                        if d4.endswith('methylation_analysis.txt'):
                            patient=d4.split('__')[2][:15]
                            values=[]
                            inf=open(path+d+'/DNA_Methylation/'+d3+'/Level_3/'+d4,'rU')
                            inf.readline()
                            for line in inf:
                                li=line.strip().split('\t')
                                if len(li)>2:
                                    values.append(li[2])
                                else:
                                    pass
                            inf.close()
                            if (patient,batchid) in methyl.keys():
                                methyl[(patient,batchid)].append(values)
                            else:
                                methyl[(patient,batchid)]=[values]
                else:
                    pass
        else:
            pass

    else:
        pass

#################Caculate the mean of the replicates!!##############################

methyl_final={}

for Agi in methyl.keys():
    if len(methyl[Agi]) >1:
        value_final=[]
        for i in range(0,len(methyl[Agi][0])):
            sum=0
            remove=0
            for j in range(0,len(methyl[Agi])):
                if methyl[Agi][j][i]=='NA':
                    remove=remove+1
    
                else:
                    sum=sum+float(methyl[Agi][j][i])
            value_final.append(sum/float(len(methyl[Agi])-remove))
        #print value_final
        methyl_final[Agi]=value_final
    else:
        methyl_final[Agi]=RNAseq[Agi][0]
methyl=0

######################The genes name list&&!!!Modify the sample path!!!!###########################    
methyl_sample='../COAD/COAD_Batch76_level3/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/jhu-usc.edu__HumanMethylation450__TCGA-A6-2671-01A-01D-1407-05__methylation_analysis.txt'
inf2=open(methyl_sample,'rU')
methyl_genes=[]
inf2.readline()
for line2 in inf2:
    li2=line2.split('\t')
    if len(li2)>2:
        methyl_genes.append('_'.join(li2[3:6]))
    else:
        pass
inf2.close()

########################Write the combat input files##########################
methyl_exp=open(path+'methylall.txt','w')
methyl_info=open(path+'methylsample.txt','w')
methyl_info.write('Array name\tSample name\tBatch\n')
Agi_keys=methyl_final.keys()
i=0
for key in Agi_keys:
    if i==len(Agi_keys)-1:
        methyl_exp.write(key[0]+'_'+key[1]+'\n')
    else:
        methyl_exp.write(key[0]+'_'+key[1]+'\t')
    i=i+1
    methyl_info.write(key[0]+'_'+key[1]+'\t'+key[0]+'_'+key[1]+'\t'+key[1]+'\n')
for i in range(0,len(methyl_genes)):
    methyl_exp.write(methyl_genes[i])
    for j in range(0,len(Agi_keys)):
        #print Agi_final
        methyl_exp.write('\t'+str(methyl_final[Agi_keys[j]][i]))
    methyl_exp.write('\n')

methyl_exp.close()
methyl_info.close()
