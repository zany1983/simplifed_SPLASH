#! /usr/bin/env python3
import argparse
import re,sys,os,math
import numpy as np
import pandas as pd

__author__='xiedejian'
__mail__='xiedejian@frasergen.com'

def read_options():
	parser=argparse.ArgumentParser(description=__doc__,
	formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-i','--insulation',help='insulation',dest='insulation',type=str,required=True)
	parser.add_argument('-nt','--noise',help='boundstrength cutoff',dest='noise',type=float,required=True)
	parser.add_argument('-o','--outdir',help='outdir',dest='outdir',type=str,required=True)
	args=parser.parse_args()
	return args

def main():
	args=read_options()
	insulation_file=os.path.abspath(args.insulation)
	if not os.path.exists(insulation_file):
		from glob import glob
		insulation_files=glob(insulation_file)
		if len(insulation_files)==0:
			print("there's no insulation files there ,please check!!!")
		elif len(insulation_files)>1:
			print('there are one more insualtion files,very confusing!!!')
		else:
			insulation_file=insulation_files[0]
	ibname=os.path.basename(insulation_file)
	noise=args.noise
	outdir=os.path.abspath(args.outdir)
	df=pd.read_table(insulation_file,header=0)
	df=df.fillna(0)
	bs=np.array(df['boundaryStrength'])
	dbs=np.diff(bs)
	splitIndex= np.argwhere(dbs!=0)
	SI=list((splitIndex+1).ravel())
	AI=[0,]
	AI.extend(SI[:-1])
	groups=list(zip(AI,SI))
	boundaries=colnames=[]
	for i,g in enumerate(groups):
		sub_group=df.iloc[g[0]:g[1],:]
		sub_group=sub_group.reset_index(drop=True)
		max=np.argmax(np.abs(sub_group['insulationScore']))
		colnames=list(sub_group.columns)
		boundaries.append(list(sub_group.iloc[max,:]))
	bounds=pd.DataFrame(boundaries,columns=colnames)
	bounds=bounds[bounds['boundaryStrength']>noise]
	bounds=bounds.reset_index(drop=True)
	out_refine_bound="{}/{}.refine_boundaries".format(outdir,ibname)
	bounds.to_csv(out_refine_bound,sep='\t',index=False)

if __name__=='__main__':
	main()
