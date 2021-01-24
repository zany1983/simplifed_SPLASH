'''
    This script is for iced_asis_matrix convert to npy cis_matrix and filterZero(for 3Dstructure)
'''

#!/usr/bin/env python
# coding: utf-8 -*- 
import argparse
import pandas as pd
import numpy as np
import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
plt.style.use('ggplot')
#import itertools
#from functools import reduce
from scipy import sparse

pd.set_option('display.precision', 2)

__author__ = 'Xiedejian'
__mail__ = 'xiedejian@frasergen.com'
__date__ = '20210123'



def main():
    parser=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
    parser.add_argument('-m','--icedAsisMatrix',help='iced asis matrix',dest='icedAsisMatrix',type=str,required=True)
    parser.add_argument('-b','--bed',help='abs_bed',dest='bed',type=str,required=True)
    parser.add_argument('-c','--chrom',help='chrom eg: chr3',dest='chrom',type=str,required=True)
    parser.add_argument('-s','--sample',help='sample',dest='sample',type=str,required=True)
    parser.add_argument('-o','--outDir',help='outDir',dest='outDir',type=str,required=True)

    args=parser.parse_args()
    
    #matrix_file, bed_file, sample, chrom, outDir = argvs
    sample = args.sample
    outDir = args.outDir
    chrom = args.chrom
    #-----------------load file------------------#
    bed = pd.read_csv(args.bed,header=None,index_col = None,names=['chrom','start','end','bNum'],sep='\t',encoding='utf-8')
    bed = bed[bed['chrom'] == chrom]
    minN = bed.iloc[0,3]
    maxN = bed.iloc[-1,3]
    resolution = bed.iloc[0,2] - bed.iloc[0,1]
    #-----------------load iced asis matrix-------#
    df = pd.read_csv(args.icedAsisMatrix, header=None,index_col = None,encoding='utf-8',sep='\t')
    df = df[df.iloc[:,0].isin(bed['bNum']) & df.iloc[:,1].isin(bed['bNum'])]

    df.loc[:,0] = df.loc[:,0] - minN
    df.loc[:,1] = df.loc[:,1] - minN

    #---------------------asis to complete matrix----------------------#    
    N = maxN - minN + 1
    counts = sparse.coo_matrix((df[2], (df[0], df[1])), shape=(N, N), dtype=float)
    df = pd.DataFrame(counts.todense(),index = bed['bNum'].tolist(),columns = bed['bNum'].tolist())
    counts = None

    #-----------------------matrix filter zeros and to npy-----------------#
    outnpy = os.path.join(outDir,'{0}_{1}_{2}_cis.npy'.format(sample,resolution,chrom))
    #heatmap(df,resolution,sample,outFile,chrom)

    df = df.loc[df.any(1),df.any(0)]
    nd = np.matrix(df)
    np.set_printoptions(precision=2)
    np.save(outnpy,nd)

    #df.to_csv(outmatrix,sep="\t",index=True,header=True,float_format='%0.2f')

    outChrSize = '{0}/budding_{1}_structure'.format(outDir,chrom)
    chrsize = (nd.shape[0] - 1) * resolution
    budding = open (outChrSize,"w")
    budding.write(str(chrsize))
    budding.close()

if __name__=="__main__":    
    #argvs = ['E:\\python\\iced_asis2matrix_heatmap\\SD_YZ_H_100000_iced.matrix','E:\\python\\iced_asis2matrix_heatmap\\SD_YZ_H_100000_abs.bed','NNJ_YZ_H','chr3','E:\\python\\iced_asis2matrix_heatmap'] 
    #argvs = sys.argv[1:]
    main()
