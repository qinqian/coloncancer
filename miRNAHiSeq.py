#!/usr/bin/env python
"""

"""
# QinQian

##
import os
import glob
os.chdir("../COAD/miRNA/miRNASeq/BCGSC__IlluminaHiSeq_miRNASeq/Level_3/")
dl = glob.glob("*mirna_quanti*.txt")


def prepr(dl):
    """
    dl: datalist input
    """
    nmlist = []
    pd = {} # patient with mirna value
    for name in dl:
        nm = name.split("__")[2][:15]
        nmlist.append(nm)

        with open(name) as mirna:
            mirna.readline()
            dall = [d.split('\t')[3] for d in mirna.readlines()]
            ID = [d.split('\t')[1] for d in mirna.readlines()]
        pd[nm] = pd.get(nm, []) + [dall]
    test = {}
    for l in nmlist:
        test[l] = test.get(l, 0) + 1

    for i in pd.keys():
        if len(pd[i]) >= 2:
        #TCGA-A6-5665-01 TCGA-A6-5661-01 has batches
            print i
    return pd, nmlist, ID

def outaver(d, nmlist):
    """
    dictionary containing all data
    """
    with open("miRNAHiseq.txt", 'w') as f:
        print >>f, '\t'.join(d.keys())
        for index in d:

    f.close()
    return

pd, nmlist = prepr(dl)
outaver(pd, nmlist)
