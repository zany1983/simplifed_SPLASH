#! /usr/bin/env python3
from glob import glob
import os,re,sys,argparse
import numpy as np
import pandas as pd

__author__='xiedejian'
__mail__  ='xiedejian@frasergen.com'

def check_files(file_path):
	path=os.path.abspath(file_path)
	file=''
	if not os.path.exists(path):
		files=glob(path)
		if len(files)==0:
			sys.stderr.write('the path:\n{}\nis not vaild\n'.format(path))
			sys.exit()
		elif len(files)>1:
			sys.stderr.write('the path:\n{}\n is  confused when there are more than one files'.format(path))
			sys.exit()
		else:
			file=files[0]
	else:
		file=path
	return file

def main():
	parser=argparse.ArgumentParser(description=__doc__,
			formatter_class=argparse.RawDescriptionHelpFormatter,
			epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-b','--bound',help='boundary',dest='boundary',required=True)
	parser.add_argument('-cs','--chromsize',help='chromsize',dest='chromsize',required=True)
	parser.add_argument('-c','--chrom',help='chrom',dest='chrom',required=True)
	parser.add_argument('-s','--sample',help='sample',dest='sample',required=True)
	parser.add_argument('-r','--resolution',help='resolution',dest='resolution',required=True)
	parser.add_argument('-o','--outdir',help='outdir',dest='outdir',required=True)
	args=parser.parse_args()
	boundary=check_files(args.boundary)
	chromsize=os.path.abspath(args.chromsize)
	size={}
	with open(chromsize,'r') as cs:
		for line in cs:
			chr,chr_size=re.split('\s+',line.rstrip())
			size[chr]=chr_size
	sample=args.sample
	outdir=os.path.abspath(args.outdir)
	chrom=args.chrom
	resolution=args.resolution
	bound_df=pd.read_table(boundary,header=0)
	start_list=list(bound_df['end'])
	start_list.insert(0, 0)
	start_df=pd.Series(start_list)
	end_list=list(bound_df['start'])
	end_list.append(size[chrom])
	end_df=pd.Series(end_list)
	tad_df=pd.concat([start_df,end_df],axis=1)
	tad_df['chr']=chrom
	tad_df.columns=['start','end','chrom']
	tad_df=tad_df[['chrom','start','end']]
	#tad_df=tad_df[(tad_df['end']-tad_df['start']>400000) & (tad_df['end']-tad_df['start']<3000000)]
	#tad_df=tad_df.reset_index(drop=True)
	tad_bound=pd.merge(tad_df,bound_df[['end','boundaryStrength']],left_on='start',right_on='end',how='left')
	tad_bound.columns=['chrom','start','end','end_y','boundaryStrength']
	tad_bound=tad_bound[['chrom','start','end','boundaryStrength']]
	tad_bound=pd.merge(tad_bound,bound_df[['start','boundaryStrength']],left_on='end',right_on='start',how='left',suffixes=['_left','_right'])
	tad_bound.columns=['chrom','start','end','boundaryStrength_left','start_right','boundaryStrength_right']
	tad_bound=tad_bound[['chrom','start','end','boundaryStrength_left','boundaryStrength_right']]
	tad_df.to_csv("{}/{}_{}_{}.filefinaldomaincalls".format(outdir,sample,resolution,chrom),sep="\t",header=False,index=False)
	tad_bound.to_csv("{}/{}_{}_{}.tad_bound.strength".format(outdir,sample,resolution,chrom),sep='\t',header=True,index=False)

if __name__=='__main__':
	main()
