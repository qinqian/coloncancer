import os
import sys
import string

#path='/mnt/Storage/home/duanxk/TCGA/GBM/'
path=sys.argv[1]
dirs=os.listdir(path)
disease=path.split('/')[-2]
Agilent={}
Agilent_genes=[]
Affy={}
Affy_genes=[]
for d in dirs:
    if d.startswith(disease+'_Batch'):
        dirs2=os.listdir(path+d)
        if 'Expression-Genes' in dirs2:
            batchid=d.split('_')[1][5:]
            dirs3=os.listdir(path+d+'/Expression-Genes')
            for d3 in dirs3:
                if d3.startswith('UNC_'):
                    #os.chdir(path+d+'/Expression-Genes/'+d3+'/Level_3')
                    dirs4=os.listdir(path+d+'/Expression-Genes/'+d3+'/Level_3')
                    for d4 in dirs4:
                        if not d4.endswith('_temp') and not d4.startswith('Batch'):
                            patient=d4.split('__')[2][:15]
                            values=[]
                            inf=open(path+d+'/Expression-Genes/'+d3+'/Level_3/'+d4,'rU')
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
                            if (patient,batchid) in Agilent.keys():
                                Agilent[(patient,batchid)].append(values)
                            else:
                                Agilent[(patient,batchid)]=[values]
                elif d3.startswith('BI'):
                    dirs4=os.listdir(path+d+'/Expression-Genes/'+d3+'/Level_3')
                    for d4 in dirs4:
                        if not d4.endswith('_temp') and not d4.startswith('Batch'):
                            patient=d4.split('__')[2][:15]
                            values=[]
                            inf=open(path+d+'/Expression-Genes/'+d3+'/Level_3/'+d4,'rU')
                            inf.readline()
                            for line in inf:
                                li=line.strip().split('\t')
                                if len(li)>2:
                                    values.append(li[2])
                                else:
                                    pass
                            inf.close()
                            if (patient,batchid) in Affy.keys():
                                Affy[(patient,batchid)].append(values)
                            else:
                                Affy[(patient,batchid)]=[values]
#################Caculate the mean of the replicates!!##############################

Agi_final={}
Affy_final={}
for Agi in Agilent.keys():
    if len(Agilent[Agi]) >1:
        value_final=[]
        for i in range(0,len(Agilent[Agi][0])):
            sum=0
            remove=0
            for j in range(0,len(Agilent[Agi])):
                if Agilent[Agi][j][i]=='null':
                    remove=remove+1
    
                else:
                    sum=sum+float(Agilent[Agi][j][i])
                
            value_final.append(sum/float(len(Agilent[Agi])-remove))
        #print value_final
        Agi_final[Agi]=value_final
    else:
        Agi_final[Agi]=Agilent[Agi][0]
for Af in Affy.keys():
    if len(Affy[Af]) >1:
        value_final=[]
        for i in range(0,len(Affy[Af][0])):
            sum=0
            remove=0
            for j in range(0,len(Affy[Af])):
                if Affy[Af][j][i]=='null':
                    remove=remove+1
                else:
                    sum=sum+float(Affy[Af][j][i])
            value_final.append(sum/float(len(Affy[Af])-remove))
        Affy_final[Af]=value_final
    else:
        Affy_final[Af]=Affy[Af][0]
Agilent=0
Affy=0
######################The genes name list&&!!!Modify the sample path!!!!###########################    
Agilent_sample='../COAD/COAD_Batch76_level3/Expression-Genes/UNC__AgilentG4502A_07_3/Level_3/unc.edu__AgilentG4502A_07_3__TCGA-A6-2671-11A-01R-1758-07__gene_expression_analysis.txt'
#Affy_sample='/mnt/Storage/home/duanxk/TCGA/GBM/GBM_Batch1_level3/Expression-Genes/BI__HT_HG-U133A/Level_3/broad.mit.edu__HT_HG-U133A__TCGA-02-0001-01C-01R-0177-01__gene_expression_analysis.txt'
inf2=open(Agilent_sample,'rU')
#inf3=open(Affy_sample,'rU')
inf2.readline()
#inf3.readline()
for line2 in inf2:
    li2=line2.split('\t')
    if len(li2)>2:
        Agilent_genes.append(li2[1])
    else:
        pass
#for line3 in inf3:
#    li3=line3.split('\t')
#    if len(li3)>2:
#        Affy_genes.append(li3[1])
#    else:
#        pass
#
inf2.close()
#inf3.close()
########################Write the combat input files##########################
Affy_exp=open(path+'Affy_all_expression.txt','w')
Affy_info=open(path+'Affy_sample_batch_info.txt','w')
Agilent_exp=open(path+'Agilent_all_expression.txt','w')
Agilent_info=open(path+'Agilent_sample_batch_info.txt','w')
Agilent_info.write('Array name\tSample name\tBatch\n')
Affy_info.write('Array name\tSample name\tBatch\n')
Agi_keys=Agi_final.keys()
i=0
for key in Agi_keys:
    if i==len(Agi_keys)-1:
        Agilent_exp.write(key[0]+'_'+key[1]+'\n')
    else:
        Agilent_exp.write(key[0]+'_'+key[1]+'\t')
    i=i+1
    Agilent_info.write(key[0]+'_'+key[1]+'\t'+key[0]+'_'+key[1]+'\t'+key[1]+'\n')
for i in range(0,len(Agilent_genes)):
    Agilent_exp.write(Agilent_genes[i])
    
    for j in range(0,len(Agi_keys)):
        #print Agi_final
        Agilent_exp.write('\t'+str(Agi_final[Agi_keys[j]][i]))
    Agilent_exp.write('\n')

Af_keys=Affy_final.keys()
i=0
for key in Af_keys:
    #print key
    if i==len(Af_keys)-1:
        Affy_exp.write(key[0]+'_'+key[1]+'\n')
    else:
        Affy_exp.write(key[0]+'_'+key[1]+'\t')
    Affy_info.write(key[0]+'_'+key[1]+'\t'+key[0]+'_'+key[1]+'\t'+key[1]+'\n')
    i=i+1
for i in range(0,len(Affy_genes)):
    Affy_exp.write(Affy_genes[i])
    for j in range(0,len(Af_keys)):
        Affy_exp.write('\t'+str(Affy_final[Af_keys[j]][i]))
    Affy_exp.write('\n')
Affy_exp.close()
Affy_info.close()
Agilent_exp.close()
Agilent_info.close()


                    
                    
                            
