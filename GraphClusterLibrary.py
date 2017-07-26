"""
Script to graph the clustermap analysis, and return clustered info
Otherwise the settings mess with the other graphs

Wren Saylor
July 12 2017

To Do:
tissue x uce
tissue x meth/c
uce x meth/c
methylation flanks
"""

import argparse
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Reduce each id/tissue x location to get a corresponding tissue/id list of associations
def listOverlap(dataframe,yItem,zItem,xItem):
	# Separate by strand
	PlusMeth = dataframe.loc[dataframe['Cytosine'] == 'C']
	MinusMeth = dataframe.loc[dataframe['Cytosine'] == 'G']
	
	# Subset just columns to use
	PlussubMeth = PlusMeth[['methLoc',yItem,zItem,xItem]]
	MinussubMeth = MinusMeth[['methLoc',yItem,zItem,xItem]]

	# Grouping a collection of a values in the xItem column
	PlusgroupMeth = PlussubMeth.join(PlussubMeth.groupby(['methLoc',yItem])[xItem].unique(),on=['methLoc',yItem],rsuffix='_r')
	MinusgroupMeth = MinussubMeth.join(MinussubMeth.groupby(['methLoc',yItem])[xItem].unique(),on=['methLoc',yItem],rsuffix='_r')

	# Just the new xItem list, and location
	PlusxLMeth = PlusgroupMeth[['methLoc','{0}_r'.format(xItem)]]
	MinusxLMeth = MinusgroupMeth[['methLoc','{0}_r'.format(xItem)]]

	print 'Collected list of {0} item associated with Location'.format(xItem)

	# Just the new xItem list, and yItem
	PlusxYMeth = PlusgroupMeth[[yItem,'{0}_r'.format(xItem)]]
	MinusxYMeth = MinusgroupMeth[[yItem,'{0}_r'.format(xItem)]]
	
	print 'Collected list of {0} item associated with {1}'.format(xItem,yItem)

	# Reindex the methLocation x ID/Tissue df
# 	PlusxYMeth.set_index('methLoc',drop=True,inplace=True)
# 	MinusxYMeth.set_index('methLoc',drop=True,inplace=True)
# 	
# 	reinPlusxYMeth = PlusxYMeth.reindex(new_index,fill_value=0)
# 	reinMinusxYMeth = MinusxYMeth.reindex(new_index,fill_value=0)
	
	
	return PlusxLMeth,MinusxLMeth,PlusxYMeth,MinusxYMeth

# Transform the Frequency, Percentage and Coverage data into graphable data frames, returning just the info for the element
def elemenetIndex(dataframe,yItem,zItem,num,uce,halfwindow,window):

	# x item is methLoc, y item is either tissue or id, z item is coverage, percentage, or frequency
	new_index = range(0,num)
	
	# Separate by strand
	PlusMeth = dataframe.loc[dataframe['Cytosine'] == 'C']
	MinusMeth = dataframe.loc[dataframe['Cytosine'] == 'G']
	
	# Subset just columns to use
	PlussubMeth = PlusMeth[['methLoc',yItem,zItem]]
	MinussubMeth = MinusMeth[['methLoc',yItem,zItem]]
	
	# Sort ascending, in order to only use the highest value with keep = last
	PlussortMeth = PlussubMeth.sort_values(['methLoc'],ascending=True)
	MinussortMeth = MinussubMeth.sort_values(['methLoc'],ascending=True)

	PlusdupMeth = PlussortMeth.drop_duplicates(['methLoc',yItem,zItem],keep='last')
	MinusdupMeth = MinussortMeth.drop_duplicates(['methLoc',yItem,zItem],keep='last')
	
	# Pivot the data frame so that each tissue/cell type is a column
	PluspivotMeth = pd.pivot_table(PlusdupMeth,index='methLoc',columns=[yItem],values=zItem,fill_value=0)
	MinuspivotMeth = pd.pivot_table(MinusdupMeth,index='methLoc',columns=[yItem],values=zItem,fill_value=0)
	
	PluspivotMeth.columns.name = None
	MinuspivotMeth.columns.name = None
	
	# Give new index, using the methLocations
	PlusindexMeth = PluspivotMeth.reindex(new_index,fill_value=0)
	MinusindexMeth = MinuspivotMeth.reindex(new_index,fill_value=0)

	# Remove the index column name
	PlusindexMeth.index.name = None
	MinusindexMeth.index.name = None
	
	# Get just the element 
	Pluselement = PlusindexMeth[(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)]
	Minuselement = MinusindexMeth[(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)]

	# Transpose the data frame for easy input into the heatamp
	PlustransMeth = Pluselement.T
	MinustransMeth = Minuselement.T
	
	PlustransMeth = PlustransMeth[PlustransMeth.columns].astype(float)
	MinustransMeth = MinustransMeth[MinustransMeth.columns].astype(float)
	
	print 'Converted {0} by {1} into data frame'.format(yItem,zItem)
	
	return PlustransMeth,MinustransMeth

# Make some graphs for fangs
def graphCluster(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine,methylationflank):

	# Parameters that all graphs will use
	fillX = range(0,(num-window))
	halfwindow = ((window/2)+1)

	# Get mean and standard deviation for AT
	ATNames = [names.index(i) for i in names if 'A' in i or 'T' in i]
	ATDataFrames = [slidingWinDF[i] for i in ATNames]
	ATconcat = pd.concat(ATDataFrames,axis=1)
	ATgroup = ATconcat.groupby(ATconcat.columns,axis=1).sum()
	ATmean = ATgroup.mean()
	ATstd = ATgroup.std()
	ATelement = ATgroup.T[(((num-uce)/2)-halfwindow):(((num-uce)/2)+uce-halfwindow)]

	# Title info
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Cluster_{0}.pdf'.format(fileName))

	# Use the row_colors to color those with similar SD?
# 	lut = dict(zip(species.unique(), "rbg"))
# 	row_colors = species.map(lut)
	heatmap0 = sns.clustermap(ATelement.T,cmap='RdPu',vmin=0,vmax=100,xticklabels=100,col_cluster=False)#,row_colors=row_colors
	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap0.ax_heatmap.set_yticks([]))
	plt.setp(heatmap0.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap0.ax_heatmap.set_ylabel('{0} UCEs'.format(len(ATelement.index)),size=8))
	plt.setp(heatmap0.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap0.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap0.ax_heatmap.set_title('Mean AT Content per UCE',size=12))
	
	sns.despine()
	pp.savefig()

	# Various combinations to plot on heatmaps, just for element
	# Frequency x (tissue, id)
	FreqPlusID,FreqMinusID = elemenetIndex(pdMeth,'id','methFreq',num,uce,halfwindow,window)
	FreqPlusTis,FreqMinusTis = elemenetIndex(pdMeth,'tissue','methFreq',num,uce,halfwindow,window)
	
	PlusLMethID,MinusLMethID,PlusYMethID,MinusYMethID = listOverlap(pdMeth,'id','methFreq','tissue')
	PlusLMethTis,MinusLMethTis,PlusYMethTis,MinusYMethTis = listOverlap(pdMeth,'tissue','methFreq','id')
	
	print PlusLMethID
	print PlusYMethID
	# zip a dictionary for location and y axis, to use as row color
	
	# Remove UCEs with out methylation within the element
	FreqPlusID = FreqPlusID[(FreqPlusID.T != 0).any()]
	FreqMinusID = FreqMinusID[(FreqMinusID.T != 0).any()]


	# Make heatmap for # methylation on pos strand (Frequency)
# 	lut = dict(zip(species.unique(), "rbg"))
# 	row_colors = species.map(lut)
	ylabels1 = FreqPlusTis.index
	heatmap1 = sns.clustermap(FreqPlusTis,cmap='RdPu',xticklabels=100,col_cluster=False,yticklabels=ylabels1)#cbar_ax=cbar5_ax,vmin=0,vmax=5
	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap1.ax_heatmap.set_ylabel('Sample',size=8))
	plt.setp(heatmap1.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap1.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap1.ax_heatmap.set_yticks(np.arange(FreqPlusTis.shape[0]) + 0.6, minor=False))
	plt.setp(heatmap1.ax_heatmap.set_yticklabels(ylabels1,minor=False,rotation=0))#,minor=False
	plt.setp(heatmap1.ax_heatmap.set_title('Methylation Frequency on Plus Strand',size=12))
	
	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
# 	lut = dict(zip(species.unique(), "rbg"))
# 	row_colors = species.map(lut)
	ylabels2 = FreqMinusTis.index
	heatmap2 = sns.clustermap(FreqMinusTis,cmap='RdPu',xticklabels=100,col_cluster=False,yticklabels=ylabels2)#cbar_ax=cbar5_ax#,vmin=0,vmax=5
	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap2.ax_heatmap.set_ylabel('Sample',size=8))
	plt.setp(heatmap2.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap2.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap2.ax_heatmap.set_yticks(np.arange(FreqMinusTis.shape[0]) + 0.6, minor=False))
	plt.setp(heatmap2.ax_heatmap.set_yticklabels(ylabels2,minor=False,rotation=0))#,minor=False
	plt.setp(heatmap2.ax_heatmap.set_title('Methylation Frequency on Minus Strand',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
# 	lut = dict(zip(tissue.unique(), "rbg"))
# 	row_colors = tissue.map(lut)
	ylabels3 = FreqPlusID.index
	heatmap3 = sns.clustermap(FreqPlusID,cmap='RdPu',xticklabels=100,col_cluster=False,yticklabels=ylabels3)#,vmin=0,vmax=5
	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap3.ax_heatmap.set_ylabel('{0} UCEs'.format(len(FreqPlusID.index)),size=8))
	plt.setp(heatmap3.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap3.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap3.ax_heatmap.set_yticks(np.arange(FreqPlusID.shape[0]) + 0.6, minor=False))
	plt.setp(heatmap3.ax_heatmap.set_yticklabels(ylabels3,minor=False,rotation=0))#,minor=False
	plt.setp(heatmap3.ax_heatmap.set_title('Methylation Frequency on Plus Strand',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on neg strand (Frequency)
# 	lut = dict(zip(tissue.unique(), "rbg"))
# 	row_colors = tissue.map(lut)
	ylabels4 = FreqMinusID.index
	heatmap4 = sns.clustermap(FreqMinusID,cmap='RdPu',xticklabels=100,col_cluster=False,yticklabels=ylabels4)#,vmin=0,vmax=5
	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap4.ax_heatmap.set_ylabel('{0} UCEs'.format(len(FreqMinusID.index)),size=8))
	plt.setp(heatmap4.ax_heatmap.set_xlabel('Position',size=8))
	plt.setp(heatmap4.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap4.ax_heatmap.set_yticks(np.arange(FreqMinusID.shape[0]) + 0.6, minor=False))
	plt.setp(heatmap4.ax_heatmap.set_yticklabels(ylabels4,minor=False,rotation=0))#,minor=False
	plt.setp(heatmap4.ax_heatmap.set_title('Methylation Frequency on Minus Strand',size=12))
# 	#http://seaborn.pydata.org/examples/timeseries_bootstrapped.html

	sns.despine()
	pp.savefig()
	pp.close()

def main(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine,methylationflank):
	graphCluster(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine,methylationflank)

if __name__ == "__main__":
	main()