#!/opt/bin/python
# Time-stamp: <2011-11-04 14:15:51 sunhf>

"""Description: Draw correlation plot for many wiggle files.

Copyright (c) 2008 Tao Liu <taoliu@jimmy.harvard.edu>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Tao Liu
@contact: taoliu@jimmy.harvard.edu
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
#import re
import logging
from optparse import OptionParser
import subprocess
from CistromeAP.taolib.CoreLib.BasicStat.Func import * 
from CistromeAP.jianlib.BwReader import BwIO

try:
    from bx.bbi.bigwig_file import BigWigFile
except:
    sys.stderr.write("Need bx-python!")
    sys.exit()

# ------------------------------------
# constants
# ------------------------------------
logging.basicConfig(level=20,
                    format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                    datefmt='%a, %d %b %Y %H:%M:%S',
                    stream=sys.stderr,
                    filemode="w"
                    )

# ------------------------------------
# Misc functions
# ------------------------------------
error   = logging.critical		# function alias
warn    = logging.warning
debug   = logging.debug
info    = logging.info
    
# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main function
# ------------------------------------
def main():
    usage = "usage: %prog <-r rfile> [options] <bigwig files> ..."
    description = "Draw correlation plot for many bigwig files. Based on qc_chIP_whole.py"
    
    optparser = OptionParser(version="%prog 0.1",description=description,usage=usage,add_help_option=False)
    optparser.add_option("-h","--help",action="help",help="Show this help message and exit.")
    #optparser.add_option("-d","--db",type="str",dest="dbname",help="UCSC db name for the assembly. Default: ce4",default="ce4")
    optparser.add_option("-r","--rfile",dest="rfile",
                         help="R output file. If not set, do not save R file.")
    optparser.add_option("-s","--step",dest="step",type="int",
                         help="sampling step in kbps. default: 100, minimal: 1",default=100)
    optparser.add_option("-z","--imgsize",dest="imgsize",type="int",
                         help="image size in inches, note the PNG dpi is 72. default: 10, minimal: 10",default=10)    
    optparser.add_option("-f","--format",dest="imgformat",type="string",
                         help="image format. PDF or PNG",default='PDF')
    #optparser.add_option("-m","--method",dest="method",type="string",default="median",
    #                     help="method to process the paired two sets of data in the sampling step. Choices are 'median', 'mean', and 'sample' (just take one point out of a data set). Default: median")
    optparser.add_option("-l","--wig-label",dest="wiglabel",type="string",action="append",
                         help="the wiggle file labels in the figure. No space is allowed. This option should be used same times as wiggle files, and please input them in the same order as -w option. default: will use the wiggle file filename as labels.")
    optparser.add_option("--min-score",dest="minscore",type="float",default=-10000,
                         help="minimum score included in calculation. Points w/ score lower than this will be discarded.")
    optparser.add_option("--max-score",dest="maxscore",type="float",default=10000,
                         help="maximum score included in calculation. Points w/ score larger than this will be discarded.")
    optparser.add_option("-H","--heatmap",dest="heatmap",action="store_true",default=False,
                         help="If True, a heatmap image will be generated instead of paired scatterplot image.")
    
    (options,wigfiles) = optparser.parse_args()

    imgfmt = options.imgformat.upper()
    if imgfmt != 'PDF' and imgfmt != 'PNG':
        print "unrecognized format: %s" % imgfmt
        sys.exit(1)

    medfunc = mean

    wigfilenum = len(wigfiles)
    if wigfilenum < 2 or not options.rfile:
        error("must provide >=2 wiggle files")
        optparser.print_help()
        sys.exit(1)

    # wig labels
    if options.wiglabel and len(options.wiglabel) == wigfilenum:
        wiglabel = options.wiglabel
    else:  # or use the filename
        wiglabel = map(lambda x:os.path.basename(x),wigfiles)
        
    if options.step < 1:
        error("Step can not be lower than 1!")
        sys.exit(1)
    if options.imgsize < 10:
        error("Image size can not be lower than 10!")
        sys.exit(1)

    # check the files
    for f in wigfiles:
        if not os.path.isfile(f):
            error("%s is not valid!" % f)
            sys.exit(1)
        
    info("number of bigwig files: %d" % wigfilenum)

    #get chromosome length from optins.wig[0]:
    p=BwIO(wigfiles[0])
    chrom_len = {}
    for i in p.chromosomeTree['nodes']:
        chrom_len[i['key']] = i['chromSize']
        
    # get the common chromosome list:
    chrset = set([t['key'] for t in p.chromosomeTree['nodes']])
    for bw in wigfiles[1:]:
        p=BwIO(bw)
        chrset = chrset.intersection(set([t['key'] for t in p.chromosomeTree['nodes']]))
    chroms = list(chrset)

    if not chroms:
        error('No common chrom found')
        sys.exit()
    info("common chromosomes are %s." % ",".join(chroms))

    # Start writing R file
    if options.rfile:
        rfhd = open(options.rfile,"w")
        rfhd.write('''require("RColorBrewer") ## from CRAN\n''')

    # for each wig file, sample...
    for i in range(len(wigfiles)):
        bw = BigWigFile(open(wigfiles[i],'rb'))
        
        info("read wiggle track from bigwig file #%d" % (i+1))
        profile = []
        for chrom in chroms:

            # The too-short chromosome will cause error in bw.summarize function below
            # So filter them out
            if chrom_len[chrom]/options.step/1000==0:
                warn("A very-short chromosome (%s) found and skipped"%chrom)
                continue
            
            summary = bw.summarize(chrom, 0, chrom_len[chrom], chrom_len[chrom]/options.step/1000)
            if not summary:
                continue
            profile_chr = summary.sum_data / summary.valid_count
            profile_chr = [str(t).replace('nan', 'NA') for t in profile_chr]
            profile.extend(profile_chr)
            
        info("write values to r file")
        rfhd.write("p%d <- c(%s)\n" %(i, ','.join(profile)))
        
    rfhd.write("c <- cbind(p0")
    for i in range(wigfilenum-1):
        rfhd.write(",p%d" % (i+1))
    rfhd.write(")\n")
    
    rfhd.write("c <- c[ c[,1]<=%f & c[,1]>=%f " % (options.maxscore,options.minscore))
    for i in range(wigfilenum-1):
        rfhd.write("& c[,%d]<=%f & c[,%d]>=%f " % (i+2,options.maxscore,i+2,options.minscore))
    rfhd.write(",]\n")
    if imgfmt == 'PDF':
        rfhd.write("pdf(\"%s.pdf\",width=%d,height=%d)\n" % (options.rfile,options.imgsize,options.imgsize))
    elif imgfmt == 'PNG':
        rfhd.write("png(\"%s.png\",units=\"in\",res=150,width=%d,height=%d)\n" % (options.rfile,options.imgsize,options.imgsize))

    if options.heatmap:                 # heatmap
        rfhd.write('library(gplots)\n')
        rfhd.write('''
m <- cor(c, method="pearson", use="pairwise.complete.obs")
''')
        labels = ",".join(map(lambda x:"\""+x+"\"",wiglabel))
        rfhd.write("rownames(m) <- c(%s)\n" % labels)
        rfhd.write("colnames(m) <- c(%s)\n" % labels)         
        rfhd.write('# draw the heatmap using gplots heatmap.2\n') 
        rfhd.write('mn <- -1\n')
        rfhd.write('mx <- 1\n')
        rfhd.write('n <- 98\n')
        rfhd.write('bias <- 1\n')
        rfhd.write('mc <- matrix(as.character(round(m, 2)), ncol=dim(m)[2])\n')
        rfhd.write('breaks <- seq(mn, mx, (mx-mn)/(n))\n')
        rfhd.write('cr <- colorRampPalette(colors = c("#2927FF","#FFFFFF","#DF5C5C"), bias=bias)\n')
        rfhd.write('heatmap.2(m, col = cr(n), breaks=breaks, trace="none", cellnote=mc, notecol="black", notecex=1.8, keysize=0.5, density.info="histogram", margins=c(27.0,27.0), cexRow=2.20, cexCol=2.20, revC=T, symm=T)\n')
    else:                               # scatterplot
        rfhd.write('''
panel.plot <- function( x,y, ... )
{
  par(new=TRUE)
  m <- cbind(x,y)
  plot(m,col=densCols(m),pch=20)
  lines(lowess(m[!is.na(m[,1])&!is.na(m[,2]),]),col="red")  
}
    
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,use="complete.obs")
  txt <- format(round(r,2),width=5,nsmall=2)
  #format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex.cor * abs(r))
  text(0.5, 0.5, txt, cex = cex.cor)
}
''')
        labels = ",".join(map(lambda x:"\""+x+"\"",wiglabel))
        rfhd.write('''
pairs(c, lower.panel=panel.plot, upper.panel=panel.cor, labels=c(%s))
''' % (labels))

    rfhd.write("dev.off()\n")
    rfhd.close()

    # try to call R
    try:
        subprocess.call(['Rscript',options.rfile])
    except:
        info("Please check %s" % options.rfile)
    else:
        info("Please check %s" % (options.rfile+'.'+imgfmt))

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
