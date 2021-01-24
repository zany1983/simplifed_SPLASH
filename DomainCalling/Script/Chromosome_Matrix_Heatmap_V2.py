#! /usr/bin/env python
import argparse
import sys
import time
import os
import copy
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.style.use('ggplot')
from scipy import sparse
import matplotlib.colors as colors
from matplotlib import collections  as mc
pd.set_option('display.precision', 2)

__author__='Xie Dejian'
__mail__= 'xiedejian@frasergen.com'
__date__= '20190929'



def colorbar_sch(axm,im,cmap,vmax):
	axins1 = inset_axes(axm, width="3%", height="6%", loc=2, bbox_to_anchor=(1.03, -0.2, 1, 1), bbox_transform=axm.transAxes,)
	ticks = [1.5]
	ticks2 = [str("= {0}".format(int(vmax)))]
	print(ticks2)
	cbar=axins1.pcolormesh([[0],[vmax]],cmap=cmap)
	axins1.set_yticks(ticks)
	axins1.set_yticklabels(ticks2,rotation=0,fontsize=12)
	axins1.set_xticks([])
	axins1.set_xticklabels([])
	axins1.yaxis.set_ticks_position('right')
	axins1.tick_params(bottom =False,top=False,left=False,right=False)


def ppm(matrix):
	df = pd.DataFrame(np.ravel(np.log(matrix)))
	df = df[df != np.inf]
	df = df[df != -np.inf]
	pp = np.nanpercentile(df, 99)
	pm = np.nanpercentile(df, 20)
	return pm, pp


def heatmap(matrix,sample,chrom,outdir,resolution,taddf):
	fig    = plt.figure(figsize=(8, 8))
	ax     = fig.add_axes([0.12,0.12,0.75,0.75])
	#cmaps = [plt.cm.get_cmap('Reds'),'YlOrBr','RdYlBu','YlOrRd','autumn_r','hot_r']
	#cmaps  = ['Greens','Blues']
	cmaps = ['Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds','YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu','GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
	pm, pp = ppm(matrix)
	matrix = matrix.replace(0,0.0001)
	#for cmap in cmaps:
	outFile   = os.path.join(outdir,'{0}_{1}_{2}_{3}_Single_Chromosome_Heatmap.png'.format(sample,resolution,chrom,cmaps[3]))
	cax    = ax.matshow(np.log(matrix),cmap=cmaps[3],clim=(pm,pp))
	ax.set_title('{0}\n{1} Cis-Matrix Heatmap'.format(sample,chrom),fontsize = 'large')
	ax.xaxis.set_ticks_position('bottom')
	ax.tick_params(bottom =True,top=False,left=True,right=False)
	ax.grid(False)
	resolution = int(resolution)
	final_list = []
	for i in taddf.index:
		tmp1_list = []
		tmp2_list = []
		start     = int(int(taddf.loc[i,'start'])/resolution)-1
		end       = int(int(taddf.loc[i,'end'])/resolution)-1
		if (start < 0):
			start=0
		tmp1_list.append((start,start))
		tmp1_list.append((end,start))
		tmp2_list.append((end,end))
		tmp2_list.append((end,start))
		final_list.append(tmp1_list)
		final_list.append(tmp2_list)
		lines = final_list
		#print lines
		lc = mc.LineCollection(lines, linewidths=0.5,color='blue')
		ax.add_collection(lc)
		ax.autoscale()
		ax.margins(0.1)
	for i in ['bottom','left','top','right']:
		ax.spines[i].set_color('black')
		ax.spines[i].set_linewidth(0.5)
		#plt.title('{0}\n{1} Cis-Matrix Heatmap'.format(sample,chrom),fontsize = 'large')
	colorbar_sch(ax,cax,cmaps[3],pp)
	fig.savefig(outFile)



def main():
	start = time.perf_counter()
	parser=argparse.ArgumentParser(description=__doc__,
	formatter_class=argparse.RawDescriptionHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}'.format(__author__,__mail__))
	parser.add_argument('-s','--sample',help='sample',dest='sample',type=str,required=True)
	parser.add_argument('-m','--matrix',help='codense Matrix',dest='matrix',type=str,required=True)
	parser.add_argument('-r','--resolution',help='resolution',dest='resolution',type=str,required=True)
	parser.add_argument('-c','--chrom',help='chrom,example:chr1',dest='chrom',type=str,required=True)
	parser.add_argument('-t','--tad',help='tadfile',dest='tad',type=str,required=True)
	parser.add_argument('-o','--outdir',help='outdir',dest='outdir',type=str,required=True)
	args  = parser.parse_args()
	#-----------------------------------Split Line---------------------------------------#
	matrix     = os.path.abspath(args.matrix)
	outdir     = os.path.abspath(args.outdir)
	tad        = os.path.abspath(args.tad)
	taddf      = pd.read_csv(tad,sep='\t',header=None,index_col=False,encoding='utf-8',names=['chrom','start','end'])
	matrixdf   = pd.read_csv(matrix,sep='\t',header=0,index_col=0,encoding='utf-8')
	sample     = args.sample
	resolution = args.resolution
	chrom      = args.chrom
	
	heatmap(matrixdf,sample,chrom,outdir,resolution,taddf)
	end     = time.perf_counter()
	runtime = end-start
	print ("this program run {0} second".format(runtime))

if __name__ == '__main__':
	main()
