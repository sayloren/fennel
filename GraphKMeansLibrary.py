"""
Script to graph k means

Wren Saylor
August 4 2017

To Do:
Doesn't work yet, need to convert list of list to array fo list
Look at other element clustering to get params, and compare (how similar are other groups to uces?)
http://scikit-learn.org/stable/modules/clustering.html#dbscan
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
from matplotlib.backends.backend_pdf import PdfPages
from GraphFangLibrary import collectDiNuc
from scipy import cluster
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

def graphKmean(dfWindow,ranWindow,names,fileName,num,uce,inuce,window,nucLine,methylationflank):

	# Parameters that all graphs will use
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)
	plt.figure(figsize=(7,7))

	# Get group, mean and standard deviation for AT
	ATgroup,ATmean,ATstd = collectDiNuc(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collectDiNuc(ranWindow,names,'A','T')
	ATelement = ATgroup.T[(((num-uce)/2)-halfwindow-methylationflank):(((num-uce)/2)+uce-halfwindow+methylationflank)]
	ranATelement = ranATgroup.T[(((num-uce)/2)-halfwindow-methylationflank):(((num-uce)/2)+uce-halfwindow+methylationflank)]
	print 'Extracted just element and methylation flank, size {0}'.format(len(ATelement))

	# Title info
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Kmeans_{0}.pdf'.format(fileName))
	sns.set_palette("husl",n_colors=8)
	
	# get the average first/last inset
	upinset = ATgroup.T[(((num-uce)/2)-halfwindow):(((num-uce)/2)-halfwindow+inuce)].mean()
	downinset = ATgroup.T[(((num-uce)/2)+uce-halfwindow-inuce):(((num-uce)/2)+uce-halfwindow)].mean()

	ATlist = [list(a) for a in zip(upinset,downinset)]
	ATarray = np.array(ATlist)
	print ATarray
	
	#https://stackoverflow.com/questions/42398403/python-k-means-clustering-array
# 	kmeans = KMeans(n_clusters=5,random_state=0).fit(ATarray)
	# see labels
	#print kmeans.labels_
	# predict new points
	#kmeans.predict([],[])
	# see where the centres of clusters are
	#kmeans.cluster_centers_
	
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	
# 	#https://stats.stackexchange.com/questions/9850/how-to-plot-data-output-of-clustering
	ATstdinitial = [cluster.vq.kmeans(ATarray,i) for i in range(1,10)]
	ax0 = plt.subplot(gs[0,:])
	ax0.plot([var for (cent,var) in ATstdinitial])
	cent,var = ATstdinitial[3]
	assignment,cdist = cluster.vq.vq(ATarray,cent)
	ax1 = plt.subplot(gs[1,:])
	ax1.scatter(ATarray[:,0],ATarray[:,1],c=assignment)

def main(dfWindow,ranWindow,names,fileName,num,uce,inuce,window,nucLine,methylationflank):
	graphKmean(dfWindow,ranWindow,names,fileName,num,uce,inuce,window,nucLine,methylationflank)

if __name__ == "__main__":
	main()