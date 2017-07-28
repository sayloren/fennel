"""
Script to graph the clustermap analysis, and return clustered info
Otherwise the settings mess with the other graphs

Wren Saylor
July 12 2017

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
import matplotlib.cbook
from scipy.spatial import distance
from scipy.cluster import hierarchy

def savePanda(pdData, strFilename):
	pdData.to_csv(strFilename, sep='\t', header=False, index=False)

# Reduce each id/tissue x location to get a corresponding tissue/id list of associations
def listOverlap(dataframe,yItem,xItem,num,uce,halfwindow,window,methylationflank):
	# check re-methFreq process

	# Separate by strand
	PlusMeth = dataframe.loc[(dataframe['Cytosine'] == 'C') & (dataframe['methLoc'] >= (((num-uce)/2)-halfwindow-methylationflank)) & (dataframe['methLoc'] <= (((num-uce)/2)+uce-halfwindow+methylationflank))]
	MinusMeth = dataframe.loc[(dataframe['Cytosine'] == 'G') & (dataframe['methLoc'] >= (((num-uce)/2)-halfwindow-methylationflank)) & (dataframe['methLoc'] <= (((num-uce)/2)+uce-halfwindow+methylationflank))]
	
	# Subset just columns to use
	PlussubMeth = PlusMeth[['methLoc',yItem,xItem]]
	MinussubMeth = MinusMeth[['methLoc',yItem,xItem]]

	PlussubMeth['methFreq'] = PlussubMeth.groupby(['methLoc',yItem,xItem])['methLoc'].transform('count')
	MinussubMeth['methFreq'] = MinussubMeth.groupby(['methLoc',yItem,xItem])['methLoc'].transform('count')

	# Grouping a collection of a values in the xItem column
	PlusgroupMeth = PlussubMeth.join(PlussubMeth.groupby(['methLoc',yItem])[xItem].unique(),on=['methLoc',yItem],rsuffix='_r')
	MinusgroupMeth = MinussubMeth.join(MinussubMeth.groupby(['methLoc',yItem])[xItem].unique(),on=['methLoc',yItem],rsuffix='_r')

	# Just the new xItem list, and location
	PlusxMeth = PlusgroupMeth[['methLoc',yItem,'{0}_r'.format(xItem)]].drop_duplicates(['methLoc',yItem],keep='last')
	MinusxMeth = MinusgroupMeth[['methLoc',yItem,'{0}_r'.format(xItem)]].drop_duplicates(['methLoc',yItem],keep='last')

	print 'Collected list of {0}s associated with Location and {1}'.format(xItem,yItem)
	
	return PlusxMeth,MinusxMeth

# Get Tissue x Id
def elementID(dataframe,yItem,zItem,num,uce,halfwindow,window,methylationflank):

	# Separate by strand
	PlusMeth = dataframe.loc[(dataframe['Cytosine'] == 'C') & (dataframe['methLoc'] >= (((num-uce)/2)-halfwindow-methylationflank)) & (dataframe['methLoc'] <= (((num-uce)/2)+uce-halfwindow+methylationflank))]
	MinusMeth = dataframe.loc[(dataframe['Cytosine'] == 'G') & (dataframe['methLoc'] >= (((num-uce)/2)-halfwindow-methylationflank)) & (dataframe['methLoc'] <= (((num-uce)/2)+uce-halfwindow+methylationflank))]
	
	# Subset just columns to use
	PlussubMeth = PlusMeth[[yItem,zItem]]
	MinussubMeth = MinusMeth[[yItem,zItem]]
	
	PlussubMeth['methFreq'] = PlussubMeth.groupby([yItem,zItem])[yItem].transform('count')
	MinussubMeth['methFreq'] = MinussubMeth.groupby([yItem,zItem])[yItem].transform('count')
	
	# Sort ascending, in order to only use the highest value with keep = last
	PlussortMeth = PlussubMeth.sort_values(['methFreq'],ascending=True)
	MinussortMeth = MinussubMeth.sort_values(['methFreq'],ascending=True)

	PlusdupMeth = PlussortMeth.drop_duplicates(['methFreq',yItem,zItem],keep='last')
	MinusdupMeth = MinussortMeth.drop_duplicates(['methFreq',yItem,zItem],keep='last')

	# Pivot the data frame so that each tissue/cell type is a column
	PluspivotMeth = pd.pivot_table(PlusdupMeth,index=[yItem],columns=[zItem],values='methFreq',fill_value=0)
	MinuspivotMeth = pd.pivot_table(MinusdupMeth,index=[yItem],columns=[zItem],values='methFreq',fill_value=0)
	
	PluspivotMeth.columns.name = None
	MinuspivotMeth.columns.name = None

	# Remove the index column name
	PluspivotMeth.index.name = None
	MinuspivotMeth.index.name = None
	
	PlusfloatMeth = PluspivotMeth[PluspivotMeth.columns].astype(float)
	MinusfloatMeth = MinuspivotMeth[MinuspivotMeth.columns].astype(float)

	print 'Collected {0} by {1} into data frame for Frequency'.format(yItem,zItem)
	
	return PlusfloatMeth,MinusfloatMeth

# Transform the Frequency, Percentage and Coverage data into graphable data frames, returning just the info for the element
def elemenetIndex(dataframe,yItem,num,uce,halfwindow,window,methylationflank):

	# x item is methLoc, y item is either tissue or id, z item is coverage, percentage, or frequency
	new_index = range(0,num)
	
	# Separate by strand
	PlusMeth = dataframe.loc[dataframe['Cytosine'] == 'C']
	MinusMeth = dataframe.loc[dataframe['Cytosine'] == 'G']
	
	# Subset just columns to use
	PlussubMeth = PlusMeth[['methLoc',yItem]]
	MinussubMeth = MinusMeth[['methLoc',yItem]]

	PlussubMeth['methFreq'] = PlussubMeth.groupby(['methLoc',yItem])['methLoc'].transform('count')
	MinussubMeth['methFreq'] = MinussubMeth.groupby(['methLoc',yItem])['methLoc'].transform('count')

	# Sort ascending, in order to only use the highest value with keep = last
	PlussortMeth = PlussubMeth.sort_values(['methLoc'],ascending=True)
	MinussortMeth = MinussubMeth.sort_values(['methLoc'],ascending=True)

	PlusdupMeth = PlussortMeth.drop_duplicates(['methLoc',yItem,'methFreq'],keep='last')
	MinusdupMeth = MinussortMeth.drop_duplicates(['methLoc',yItem,'methFreq'],keep='last')
	
	# Pivot the data frame so that each tissue/cell type is a column
	PluspivotMeth = pd.pivot_table(PlusdupMeth,index='methLoc',columns=[yItem],values='methFreq',fill_value=0)
	MinuspivotMeth = pd.pivot_table(MinusdupMeth,index='methLoc',columns=[yItem],values='methFreq',fill_value=0)
	
	PluspivotMeth.columns.name = None
	MinuspivotMeth.columns.name = None
	
	# Give new index, using the methLocations
	PlusindexMeth = PluspivotMeth.reindex(new_index,fill_value=0)
	MinusindexMeth = MinuspivotMeth.reindex(new_index,fill_value=0)

	# Remove the index column name
	PlusindexMeth.index.name = None
	MinusindexMeth.index.name = None
	
	# Get just the element 
	Pluselement = PlusindexMeth[(((num-uce)/2)-halfwindow-methylationflank):(((num-uce)/2)+uce-halfwindow+methylationflank)]
	Minuselement = MinusindexMeth[(((num-uce)/2)-halfwindow-methylationflank):(((num-uce)/2)+uce-halfwindow+methylationflank)]

	# Transpose the data frame for easy input into the heatamp
	PlustransMeth = Pluselement.T
	MinustransMeth = Minuselement.T
	
	PlustransMeth = PlustransMeth[PlustransMeth.columns].astype(float)
	MinustransMeth = MinustransMeth[MinustransMeth.columns].astype(float)
	
	print 'Converted {0} by Frequency into data frame'.format(yItem)
	
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
	ATelement = ATgroup.T[(((num-uce)/2)-halfwindow-methylationflank):(((num-uce)/2)+uce-halfwindow+methylationflank)]
	ATstd = ATelement.std()
	ATQcut = pd.qcut(ATstd,q=8,labels=False)

	# Title info
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Cluster_{0}.pdf'.format(fileName))

	# Use the row_colors to color those with similar SD?
	huslPalette = sns.husl_palette(8, s=.45)
# 	lut = dict(zip(ATQcut.unique(), huslPalette))
# 	col_colors = ATQcut.map(lut)
	heatmap0 = sns.clustermap(ATelement.T,cmap='RdPu',vmin=0,vmax=100,xticklabels=50,col_cluster=False)#,row_colors=row_colors
	plt.setp(heatmap0.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap0.ax_heatmap.set_yticks([]))
	plt.setp(heatmap0.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap0.ax_heatmap.set_ylabel('{0} UCEs'.format(len(ATelement.T.index)),size=8))
	plt.setp(heatmap0.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap0.ax_heatmap.tick_params(labelsize=10))
# 	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
# 	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap0.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap0.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.set_title('Mean AT Content per UCE',size=12))
	ATOrdered = heatmap0.dendrogram_row.linkage.reordered_ind()
	
	sns.despine()
	pp.savefig()

	# Various combinations to plot on heatmaps, just for element plus methylation flanks
	# Frequency x Tissue x ID X Location
	FreqPlusID,FreqMinusID = elemenetIndex(pdMeth,'id',num,uce,halfwindow,window,methylationflank)
	FreqPlusTis,FreqMinusTis = elemenetIndex(pdMeth,'tissue',num,uce,halfwindow,window,methylationflank)
	XPlus,XMinus = elementID(pdMeth,'id','tissue',num,uce,halfwindow,window,methylationflank)

# 	# Get the aggregated list of ids/tissues corresponding to location and tissues/ids
# 	PlusMethID,MinusMethID = listOverlap(pdMeth,'id','tissue',num,uce,halfwindow,window,methylationflank)
# 	PlusMethTis,MinusMethTis = listOverlap(pdMeth,'tissue','id',num,uce,halfwindow,window,methylationflank)

	# Remove UCEs with out methylation within the element - only for ID group
	FreqPlusID = FreqPlusID[(FreqPlusID.T != 0).any()]
	FreqMinusID = FreqMinusID[(FreqMinusID.T != 0).any()]

	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap1 = sns.clustermap(FreqPlusTis,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels1 = heatmap1.ax_heatmap.get_yticklabels()
	plt.setp(heatmap1.ax_heatmap.set_yticklabels(ylabels1,rotation=0))
	plt.setp(heatmap1.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap1.ax_heatmap.set_ylabel('Sample',size=10))
	plt.setp(heatmap1.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap1.ax_heatmap.tick_params(labelsize=10))
# 	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
# 	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap1.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap1.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.set_title('Methylation Frequency on Plus Strand',size=12))
	
	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap2 = sns.clustermap(FreqMinusTis,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels2 = heatmap2.ax_heatmap.get_yticklabels()
	plt.setp(heatmap2.ax_heatmap.set_yticklabels(ylabels2,rotation=0))
	plt.setp(heatmap2.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap2.ax_heatmap.set_ylabel('Sample',size=10))
	plt.setp(heatmap2.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap2.ax_heatmap.tick_params(labelsize=10))
# 	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
# 	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap2.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap2.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.set_title('Methylation Frequency on Minus Strand',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap3 = sns.clustermap(FreqPlusID,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels3 = heatmap3.ax_heatmap.get_yticklabels()
	plt.setp(heatmap3.ax_heatmap.set_yticklabels(ylabels3,rotation=0))
	plt.setp(heatmap3.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap3.ax_heatmap.set_ylabel('{0} UCEs'.format(len(FreqPlusID.index)),size=10))
	plt.setp(heatmap3.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap3.ax_heatmap.tick_params(labelsize=10))
# 	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
# 	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap3.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap3.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.set_title('Methylation Frequency on Plus Strand',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on neg strand (Frequency)
	heatmap4 = sns.clustermap(FreqMinusID,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels4 = heatmap4.ax_heatmap.get_yticklabels()
	plt.setp(heatmap4.ax_heatmap.set_yticklabels(ylabels4,rotation=0))
	plt.setp(heatmap4.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap4.ax_heatmap.set_ylabel('{0} UCEs'.format(len(FreqMinusID.index)),size=10))
	plt.setp(heatmap4.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap4.ax_heatmap.tick_params(labelsize=10))
# 	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+uce),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
# 	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+inuce),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap4.ax_heatmap.axvline(x=(((num-uce)/2)+(uce-inuce)),linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
# 	plt.setp(heatmap4.ax_heatmap.axvline(x=((num-uce)/2),linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.set_title('Methylation Frequency on Minus Strand',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap5 = sns.clustermap(XPlus,cmap='RdPu')
	ylabels5 = heatmap5.ax_heatmap.get_yticklabels()
	plt.setp(heatmap5.ax_heatmap.set_yticklabels(ylabels5,rotation=0))
	plt.setp(heatmap5.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap5.ax_heatmap.set_ylabel('{0} UCEs'.format(len(FreqMinusID.index)),size=10))
	plt.setp(heatmap5.ax_heatmap.set_xlabel('Sample',size=10))
	plt.setp(heatmap5.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap5.ax_heatmap.set_title('Methylation Frequency on Plus Strand',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on neg strand (Frequency)
	heatmap6 = sns.clustermap(XMinus,cmap='RdPu')
	ylabels6 = heatmap6.ax_heatmap.get_yticklabels()
	plt.setp(heatmap6.ax_heatmap.set_yticklabels(ylabels6,rotation=0))
	plt.setp(heatmap6.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap6.ax_heatmap.set_ylabel('{0} UCEs'.format(len(FreqMinusID.index)),size=10))
	plt.setp(heatmap6.ax_heatmap.set_xlabel('Sample',size=10))
	plt.setp(heatmap6.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap6.ax_heatmap.set_title('Methylation Frequency on Minus Strand',size=12))
	#heatmap6.dendrogram_col.linkage.reordered_ind
	#heatmap6.dendrogram_row.linkage.reordered_ind

	sns.despine()
	pp.savefig()
	pp.close()
	
	return ATOrdered

def main(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine,methylationflank):
	ATOrdered = graphCluster(slidingWinDF,pdMeth,names,fileName,num,uce,inuce,window,nucLine,methylationflank)
	return ATOrdered

if __name__ == "__main__":
	main()