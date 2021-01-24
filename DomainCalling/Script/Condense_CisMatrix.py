#! /usr/bin/env python
import argparse
import sys
import time
import os
import numpy as np
import pandas as pd
from scipy import sparse

__author__= 'Xie Dejian'
__mail__  = 'xiedejian@frasergen.com'
__date__  = '20190927'

def main():
	start           = time.perf_counter()
	parser          = argparse.ArgumentParser(description=__doc__,
	formatter_class = argparse.RawDescriptionHelpFormatter,
	epilog          = 'author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-s','--sample',help='sample',dest='sample',type=str,required=True)
	parser.add_argument('-b','--bedFrag',help='abs bed',dest='bedFrag',type=str,required=True)
	parser.add_argument('-m','--matrix',help='cisSparse Matrix',dest='matrix',type=str,required=True)
	parser.add_argument('-g','--genome',help='hg19',dest='genome',type=str,required=True)
	parser.add_argument('-c','--chr',help='chrom',dest='chrom',type=str,required=True)
	parser.add_argument('-o','--outdir',help='outdir',dest='outdir',type=str,required=True)
	parser.add_argument('-asym','--asymmetric',help='if need symmetric yes else no',dest='asymmetric',default="no",type=str)
	args       = parser.parse_args()
	#-------------------------------------Matrix-------------------------------------#
	bedfile    = os.path.abspath(args.bedFrag)
	matrixfile = os.path.abspath(args.matrix)
	outdir     = os.path.abspath(args.outdir)
	beddf      = pd.read_csv(bedfile,header=None,index_col=False,sep='\t',names=['chrom','start','end','bnum'],encoding='utf-8')
	df         = pd.read_csv(matrixfile,header=None,index_col=False,sep='\t',names=['b1','b2','contact'],encoding='utf-8')
	beddf      = beddf[beddf['chrom']==args.chrom]
	beddf      = beddf.drop_duplicates()
	beddf      = beddf.reset_index(drop=True)
	print(beddf)
	minN       = beddf.iloc[0,3]
	print(minN)
	maxN       = beddf.iloc[-1,3]
	df         = df[df.iloc[:,0] >= minN]
	df         = df[df.iloc[:,0] < maxN]
	df         = df[df.iloc[:,1] >= minN]
	df         = df[df.iloc[:,1] < maxN]
	df['b1'] = df['b1'] - minN
	df['b2'] = df['b2'] - minN
	N        = maxN-minN +1
	#print(N)
	#print(df['b1'])
	matrix_index = np.array(np.arange(df.shape[0]))
	df.index     = matrix_index
	#----------------------- tansform matrix ---------------------------------------------------#
	sparse_matrix= sparse.coo_matrix((df['contact'], (df['b1'], df['b2'])), shape=(N,N), dtype=float).toarray()
	#----------------------- build  symmetric matrix --------------------------------------------#
	if args.asymmetric == "yes":
		diag_matrix   = np.diag(np.diag(sparse_matrix))
		sparse_matrix = sparse_matrix.T + sparse_matrix 
		sparse_matrix = sparse_matrix-diag_matrix
	#----------------------- change colname rowname -----------------------------------------------#
	sparse_matrix = pd.DataFrame(sparse_matrix)
	sparse_matrix = sparse_matrix.fillna(0)
	inddf         = np.arange(N)
	headers_ref   = [args.genome for x in inddf]
	bin_num_df    = pd.Series(beddf['bnum']).apply(lambda x : str(x))
	headers_ref   = pd.Series(headers_ref)
	chromdf       = pd.Series(beddf['chrom'])
	startdf       = pd.Series(beddf['start']).apply(lambda x : str(x))
	#print(startdf.shape[0])
	#print(headers_ref)
	enddf         = pd.Series(beddf['end']).apply(lambda x : str(x))
	headers_suf   = chromdf.str.cat(startdf,sep=':')
	headers_suf   = headers_suf.str.cat(enddf,sep="-")
	headers       = bin_num_df.str.cat([headers_ref,headers_suf],sep="|")
	headers       = list(headers)
	sparse_matrix.columns=headers
	sparse_matrix.index=headers
	prefix="{0}/{1}.{2}_{3}.{2}_{3}.matrix".format(args.outdir,args.sample,args.genome,args.chrom)
	sparse_matrix.to_csv('{0}'.format(prefix),sep="\t",index=True,header=True,float_format='%0.3f')
	end     = time.perf_counter()
	runtime = end-start
	print ("This program run {0} second".format(runtime))

if __name__ == '__main__':
	main()
