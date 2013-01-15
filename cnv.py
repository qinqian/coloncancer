#!/usr/bin/env python
"""
Time-stamp: < modified by qinq :2013-01-13 21:34:18 >
info:
process TCGA tumor snp6.0 copy number segment data
1. intersect to get snp6.0 gene to copy number information
2. batch average
batch example                             patient      diff:plate
broad.mit.edu__Genome_Wide_SNP_6__TCGA-AA-3532-01A-01D-0819-01__snp_analysis.hg19.seg.txt
broad.mit.edu__Genome_Wide_SNP_6__TCGA-AA-3532-01A-01D-1549-01__snp_analysis.hg19.seg.txt
"""
# QinQian
##

import pandas
import glob
import os
import sys

from collections import defaultdict
from multiprocessing import Pool, Process, Queue

def locate(d = ".", ref_version = "hg19"):
    """ find data
    ## copy is the processed data from shell"""
    return glob.glob(os.path.join(d,"*%s*copy*" % ref_version))

def find_batch(dl):
    """ find patient in different
    batches by set function """
    patients_batch = [ i.split("__")[2][:15] for i in dl ]
    batch_number = {}
    for p in patients_batch:
        batch_number[p] = batch_number.get(p, 0) + 1
    patients_batches = [ i for i in batch_number if batch_number[i] > 1 ]
    simple_patient = [i for i in batch_number if batch_number[i] == 1]

    return simple_patient, patients_batches

def refgene_dict(ref):
    """ set refgene to 0 copy number """
    refgene = {}
    with open(ref, 'rU') as f:
        for line in f:
            gene = line.strip().split()[3]
            refgene[gene] = 0
    f.close()
    return refgene

def batch_average(d, p, ref_version, refgene, output):
    """ patients with batch data"""

    queue = Queue()

    processed = {}
    handle = lambda f: [i.strip().split() for i in open(f, 'rU').xreadlines()]
    bi = 0 ## assumed batch id
#    for p in batches:
    ## copy is the processed data from shell
    data = glob.glob(os.path.join(d, "*{0}*{1}*copy*".format(p, ref_version)))
    single_data = {} ## temp
    batch_data = map(handle, data)
    for batch in batch_data:
        for gene in batch:
            single_data[gene[1]] = single_data.get(gene[1], []) + [float(gene[0])]
    for k, v in single_data.iteritems():
        ## average batch
        single_data[k] = float(sum(v)) / len(v)
    processed[p] = single_data
    result = match(processed, refgene, p, output)
    return processed

def match(processed, refgene, p, output):
    """ refgene copy number matching"""
    result = {}
#    for p in processed.keys(): ## patients
    result[p] = {}
    for gene in refgene.keys():
        if gene in processed[p].keys(): ## genes' copy number dict
            result[p][gene] = processed[p][gene]
        else:
            result[p][gene] = 0
    ## write out

    write_tcga(result, refgene, output, p)

    return result

def write_tcga(td, tgene, output, p):
    """ write output patient, genes, 
    datatype batch average value tables
    td : treated data, average or median with NA removal or smooth
    tgene: gene information from different platform
    """
    outvalue = open(os.path.join(output, "snp_cnv%s.txt" % p), 'w')
    i=0 ## assumed batch ids
    for key in td.keys():
        if i==len(td.keys())-1:
            outvalue.write(key + '_' + str(i+1) +'\n')
        else:
            outvalue.write(key + '_' + str(i+1) +'\t')
        i += 1

    genes = tgene.keys()
    for i in range(0,len(genes)):
        outvalue.write(genes[i])
        for j in range(0,len(td.keys())):
            outvalue.write('\t'+str(td[td.keys()[j]][genes[i]]))
        outvalue.write('\n')
    outvalue.close()
        
def main():
    print "input: seg directory , reference version, output directory"
    dir = sys.argv[1]
    ref = sys.argv[2]
    output = sys.argv[3]
    data = locate(dir, ref)
    simple, multiple = find_batch(data)
    refgene = refgene_dict("./processed/hg19.bed")

    pool = Pool(processes=8) # wait
    all = simple + multiple
    for i in all:
        result = pool.apply_async(batch_average, (dir, i, ref, refgene, output))
    pool.close()
    pool.join()

    # paste to get the whole table
    # exec("result = %s" %f.read().strip())
    if result.successful():
        print result.keys()[0]
        print "succeed!"

if __name__ == "__main__":
    main()
