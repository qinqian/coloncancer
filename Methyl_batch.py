def write_tcga(td, tgene, outvalue, outinfo):
    """ write output patient, genes, 
    datatype batch average value tables
    td : treated data, average or median with NA removal or smooth
    tgene: gene information from different platform
    outvalue: output table with value
    outinfo: output information table
    """
    i=0
    for key in td.keys():
        if i==len(td.keys())-1:
            outvalue.write(key[0]+'_'+key[1]+'\n')
        else:
            outvalue.write(key[0]+'_'+key[1]+'\t')
        i=i+1
        outinfo.write(key[0]+'_'+key[1]+'\t'+key[0]+'_'+key[1]+'\t'+key[1]+'\n')

    for i in range(0,len(tgene)):
        outvalue.write(tgene[i])
        for j in range(0,len(td.keys())):
            outvalue.write('\t'+str(td[td.keys()[j]][i]))
        outvalue.write('\n')

def main(path):
    import os
    dirs=os.listdir(path)
    disease=path.split('/')[-2]
    methl27 = {}
    methl450 = {}
    for d in dirs:
        if d.startswith(disease+'_Batch'):
            dirs2=os.listdir(path+d)
            if 'DNA_Methylation' in dirs2:
                batchid=d.split('_')[1][5:]
                dirs3=os.listdir(path+d+'/DNA_Methylation')
                for d3 in dirs3:
                    if d3.startswith('JHU_USC__HumanMethylation27'):
                        dirs4=os.listdir(path+d+'/DNA_Methylation/'+d3+'/Level_3')
                        for d4 in dirs4:
                            if not d4.endswith('_temp') and not d4.startswith('Batch'):
                                patient=d4.split('__')[2][:15]
                                values=[]
                                inf=open(path+d+'/DNA_Methylation/'+d3+'/Level_3/'+d4,'rU')
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
                                if (patient,batchid) in methl27.keys():
                                    methl27[(patient,batchid)].append(values)
                                else:
                                    methl27[(patient,batchid)]=[values]
                    elif d3.startswith('JHU_USC__HumanMethylation450'):
                        dirs4=os.listdir(path+d+'/DNA_Methylation/'+d3+'/Level_3')
                        for d4 in dirs4:
                            if not d4.endswith('_temp') and not d4.startswith('Batch'):
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
                                if (patient,batchid) in methl450.keys():
                                    methl450[(patient,batchid)].append(values)
                                else:
                                    methl450[(patient,batchid)]=[values]
    #################Caculate the mean of the replicates!!##############################

    methl450_final={}
    methl27_final={}
    for m450 in methl450.keys():
        if len(methl450[m450]) >1:
            value_final=[]
            for i in range(0,len(methl450[m450][0])):
                sum = 0
                num = 0
                for j in range(0,len(methl450[m450])):
                    if methl450[m450][j][i]=='NA':
                        sum = sum + 0
                    else:
                        num += 1
                        sum=sum+float(methl450[m450][j][i])
                if num == 0:
                    value_final.append(0)
                else:
                    value_final.append(sum/float(num))
            methl450_final[m450]=value_final
        else:
            methl450_final[m450]=methl450[m450][0]
    for m27 in methl27.keys():
        if len(methl27[m27]) >1:
            value_final=[]
            for i in range(0,len(methl27[m27][0])):
                sum=0
                num=0
                for j in range(0,len(methl27[m27])):
                    if methl27[m27][j][i]=='NA':
                        sum = sum + 0
                    else:
                        num += 1
                        sum=sum+float(methl27[m27][j][i])
                if num == 0:
                    value_final.append(0)
                else:
                    value_final.append(sum/float(num))
            methl27_final[m27] = value_final
        else:
            methl27_final[m27]= methl27[m27][0]
    print "Sample Processed!"
    ######################The genes name list&&!!!Modify the sample path!!!!###########################    
    methl27_genes = []
    methl450_genes = []
    methl450_sample='/mnt/Storage/home/qinq/TCGA/COAD/COAD_Batch116_level3/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/jhu-usc.edu__HumanMethylation450__TCGA-A6-2675-01A-02D-1721-05__methylation_analysis.txt'
    methl27_sample='/mnt/Storage/home/qinq/TCGA/COAD/COAD_Batch41_level3/DNA_Methylation/JHU_USC__HumanMethylation27/Level_3/jhu-usc.edu__HumanMethylation27__TCGA-A6-3807-01A-01D-1020-05__methylation_analysis.txt'
    inf2=open(methl450_sample,'rU')
    inf3=open(methl27_sample,'rU')
    inf2.readline()
    inf3.readline()
    for line2 in inf2:
        li2=line2.split('\t')
        if len(li2)>2:
            methl450_genes.append('_'.join(li2[3:6]).strip("\n"))
    for line3 in inf3:
        li3=line3.split('\t')
        if len(li3)>2:
            methl27_genes.append('_'.join(li3[3:6]).strip("\n"))
    inf2.close()
    inf3.close()

    ## output 
    methl27_exp=open(path+'methl27_all_expression.txt','w')
    methl27_info=open(path+'methl27_sample_batch_info.txt','w')
    methl450_exp=open(path+'methl450_all_expression.txt','w')
    methl450_info=open(path+'methl450_sample_batch_info.txt','w')
    methl450_info.write('Methl450 name\tSample name\tBatch\n')
    methl27_info.write('Methl27 name\tSample name\tBatch\n')

    write_tcga(methl450_final, methl450_genes, methl450_exp, methl450_info)
    write_tcga(methl27_final, methl27_genes, methl27_exp, methl27_info)

if __name__ == "__main__": 
    import sys
    try:
        main(sys.argv[1])
    except KeyboardInterrupt:
        print >> sys.stderr, "Bye"
