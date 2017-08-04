"""
Script to graph k means

Wren Saylor
August 4 2017
"""

import argparse
import numpy as np
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from GraphFangLibrary import collectAT

def graphKmean(dfWindow,ranWindow,names,fileName,num,uce,inuce,window,nucLine,methylationflank):
# 	#https://stats.stackexchange.com/questions/9850/how-to-plot-data-output-of-clustering

	# Parameters that all graphs will use
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)
	plt.figure(figsize=(7,7))

	# Get group, mean and standard deviation for AT
	ATgroup,ATmean,ATstd = collectAT(dfWindow,names)
	ranATgroup,ranATmean,ranATstd = collectAT(ranWindow,names)
	ATelement = ATgroup.T[(((num-uce)/2)-halfwindow-methylationflank):(((num-uce)/2)+uce-halfwindow+methylationflank)]
	ranATelement = ranATgroup.T[(((num-uce)/2)-halfwindow-methylationflank):(((num-uce)/2)+uce-halfwindow+methylationflank)]
	print 'Extracted just element and methylation flank, size {0}'.format(len(ATelement))

	# Title info
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Cluster_{0}.pdf'.format(fileName))
	sns.set_palette("husl",n_colors=8)#(len(nucLine)*2)

	ATstdint = ATstd.astype(int)
	ATstdvalues = ATstdint.tolist()
	ATstdinitial = [cluster.vq.kmeans(ATstdvalues,i) for i in range(1,10)]
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	ax0 = plt.subplot(gs[0,:])
	ax0.plot([var for (cent,var) in ATstdinitial])
	cent,var = ATstdinitial[3]
	assignmetn,cdist = cluster.vq.vq(ATstdvalues,cent)
	ax1 = plt.subplot(gs[1,:])
	ax1.scatter(ATstdvalues[:,0],ATstdvalues[:,1],c=assignment)

def main(dfWindow,ranWindow,names,fileName,num,uce,inuce,window,nucLine,methylationflank):
	graphKmean(dfWindow,ranWindow,names,fileName,num,uce,inuce,window,nucLine,methylationflank)

if __name__ == "__main__":
	main()