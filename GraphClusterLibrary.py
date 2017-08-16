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
import scipy
from scipy import cluster
from scipy.spatial import distance
from scipy.cluster import hierarchy
from GraphFangLibrary import collectDiNuc
import GraphTableLibrary
import GlobalVariables

# Reduce each id/tissue x location to get a corresponding tissue/id list of associations
def listOverlap(dataframe,yItem,xItem):
	# check re-methFreq process

	# Separate by strand
	PlusMeth = dataframe.loc[(dataframe['Cytosine'] == 'C') & (dataframe['methLoc'] >= (GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank)) & (dataframe['methLoc'] <= (GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank))]
	MinusMeth = dataframe.loc[(dataframe['Cytosine'] == 'G') & (dataframe['methLoc'] >= (GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank)) & (dataframe['methLoc'] <= (GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank))]
	
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
def elementID(dataframe,yItem,zItem):

	# Separate by strand
	PlusMeth = dataframe.loc[(dataframe['Cytosine'] == 'C') & (dataframe['methLoc'] >= (GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank)) & (dataframe['methLoc'] <= (GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank))]
	MinusMeth = dataframe.loc[(dataframe['Cytosine'] == 'G') & (dataframe['methLoc'] >= (GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank)) & (dataframe['methLoc'] <= (GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank))]
	
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
def elemenetIndex(dataframe,yItem):

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
	Pluselement = PlusindexMeth[(GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank):(GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank)]
	Minuselement = MinusindexMeth[(GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank):(GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank)]

	# Transpose the data frame for easy input into the heatamp
	PlustransMeth = Pluselement.T
	MinustransMeth = Minuselement.T
	
	PlustransMeth = PlustransMeth[PlustransMeth.columns].astype(float)
	MinustransMeth = MinustransMeth[MinustransMeth.columns].astype(float)
	
	print 'Converted {0} by Frequency into data frame'.format(yItem)
	
	return PlustransMeth,MinustransMeth

# Make dictionary for row and column colors based on standard deviation
def dictColors(ATelement,huslPalette):
	ATQcutPosition = pd.qcut(ATelement.std(axis=1),q=8,labels=False)
	ATQcutElement = pd.qcut(ATelement.std(),q=8,labels=False)
	lutElement = dict(zip(ATQcutElement.unique(), huslPalette))
	elementColors = ATQcutElement.map(lutElement)
	lutPosition = dict(zip(ATQcutPosition.unique(), huslPalette))
	positionColors = ATQcutPosition.map(lutPosition)
	print 'Made dictionary for standard deviation'
	return elementColors,positionColors

# Make some graphs for fangs
def graphCluster(dfWindow,ranWindow,pdMeth,rnMeth,names,fileName):

	plt.figure(figsize=(7,7))

	# Get group, mean and standard deviation for AT
	ATgroup,ATmean,ATstd = collectDiNuc(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collectDiNuc(ranWindow,names,'A','T')
	ATelement = ATgroup.T[(GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank):(GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank)]
	ranATelement = ranATgroup.T[(GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank):(GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank)]
	print 'Extracted just element and methylation flank, size {0}'.format(len(ATelement))

	# Title info
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Cluster_{0}.pdf'.format(fileName))
	sns.set_palette("husl",n_colors=8)#(len(nucLine)*2)

# 	Use the row_colors to color those with similar SD?
	huslPalette = sns.husl_palette(8, s=.45)
	elementColors,positionColors = dictColors(ATelement,huslPalette)
	heatmap0 = sns.clustermap(ATelement.T,cmap='RdPu',vmin=0,vmax=100,xticklabels=50,col_cluster=False,row_colors=elementColors,col_colors=positionColors)
	plt.setp(heatmap0.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap0.ax_heatmap.set_yticks([]))
	plt.setp(heatmap0.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap0.ax_heatmap.set_ylabel('{0} UCEs'.format(len(ATelement.T.index)),size=8))
	plt.setp(heatmap0.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap0.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap0.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap0.ax_heatmap.set_title('Mean AT Content per Element',size=12))
# 	ATOrdered = heatmap0.dendrogram_row.reordered_ind
	
	sns.despine()
	pp.savefig()
	
# 	Use the row_colors to color those with similar SD?
	ranelementColors,ranpositionColors = dictColors(ranATelement,huslPalette)
	heatmap1 = sns.clustermap(ranATelement.T,cmap='RdPu',vmin=0,vmax=100,xticklabels=50,col_cluster=False,row_colors=ranelementColors,col_colors=ranpositionColors)
	plt.setp(heatmap1.ax_heatmap.tick_params(labelsize=8))
	plt.setp(heatmap1.ax_heatmap.set_yticks([]))
	plt.setp(heatmap1.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap1.ax_heatmap.set_ylabel('{0} UCEs'.format(len(ranATelement.T.index)),size=8))
	plt.setp(heatmap1.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap1.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap1.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap1.ax_heatmap.set_title('Mean AT Content per Random Region',size=12))
# 	ranATOrdered = heatmap1.dendrogram_row.reordered_ind

	sns.despine()
	pp.savefig()
	print 'Plotted cluster plot for mean AT content for all elements and random regions'

	# Various combinations to plot on heatmaps, just for element plus methylation flanks
	# Frequency x Tissue x ID X Location
	FreqPlusID,FreqMinusID = elemenetIndex(pdMeth,'id')
	FreqPlusTis,FreqMinusTis = elemenetIndex(pdMeth,'tissue')
	XPlus,XMinus = elementID(pdMeth,'id','tissue')

	ranFreqPlusID,ranFreqMinusID = elemenetIndex(rnMeth,'id')
	ranFreqPlusTis,ranFreqMinusTis = elemenetIndex(rnMeth,'tissue')
	ranXPlus,ranXMinus = elementID(rnMeth,'id','tissue')

	# Remove UCEs with out methylation within the element - only for ID group
	FreqPlusID = FreqPlusID[(FreqPlusID.T != 0).any()]
	FreqMinusID = FreqMinusID[(FreqMinusID.T != 0).any()]

	ranFreqPlusID = ranFreqPlusID[(ranFreqPlusID.T != 0).any()]
	ranFreqMinusID = ranFreqMinusID[(ranFreqMinusID.T != 0).any()]

	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap2 = sns.clustermap(FreqPlusTis,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels2 = heatmap2.ax_heatmap.get_yticklabels()
	plt.setp(heatmap2.ax_heatmap.set_yticklabels(ylabels2,rotation=0))
	plt.setp(heatmap2.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap2.ax_heatmap.set_ylabel('Sample',size=10))
	plt.setp(heatmap2.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap2.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap2.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap2.ax_heatmap.set_title('Methylation Frequency on Plus Strand for Elements',size=12))
	
	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap3 = sns.clustermap(FreqMinusTis,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels3 = heatmap3.ax_heatmap.get_yticklabels()
	plt.setp(heatmap3.ax_heatmap.set_yticklabels(ylabels3,rotation=0))
	plt.setp(heatmap3.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap3.ax_heatmap.set_ylabel('Sample',size=10))
	plt.setp(heatmap3.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap3.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap3.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap3.ax_heatmap.set_title('Methylation Frequency on Minus Strand for Elements',size=12))

	sns.despine()
	pp.savefig()
	print 'Plotted methylation frequency for tissue types x position, for element'

	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap4 = sns.clustermap(ranFreqPlusTis,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels4 = heatmap4.ax_heatmap.get_yticklabels()
	plt.setp(heatmap4.ax_heatmap.set_yticklabels(ylabels4,rotation=0))
	plt.setp(heatmap4.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap4.ax_heatmap.set_ylabel('Sample',size=10))
	plt.setp(heatmap4.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap4.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap4.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap4.ax_heatmap.set_title('Methylation Frequency on Plus Strand for Random Regions',size=12))
	
	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap5 = sns.clustermap(ranFreqMinusTis,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels5 = heatmap5.ax_heatmap.get_yticklabels()
	plt.setp(heatmap5.ax_heatmap.set_yticklabels(ylabels5,rotation=0))
	plt.setp(heatmap5.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap5.ax_heatmap.set_ylabel('Sample',size=10))
	plt.setp(heatmap5.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap5.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap5.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap5.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap5.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap5.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap5.ax_heatmap.set_title('Methylation Frequency on Minus Strand for Random Regions',size=12))

	sns.despine()
	pp.savefig()
	print 'Plotted methylation frequency for tissue types x position, for  random regions'

	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap6 = sns.clustermap(FreqPlusID,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels6 = heatmap6.ax_heatmap.get_yticklabels()
	plt.setp(heatmap6.ax_heatmap.set_yticklabels(ylabels6,rotation=0))
	plt.setp(heatmap6.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap6.ax_heatmap.set_ylabel('{0} Elements'.format(len(FreqPlusID.index)),size=10))
	plt.setp(heatmap6.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap6.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap6.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap6.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap6.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap6.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap6.ax_heatmap.set_title('Methylation Frequency on Plus Strand for Elements',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on neg strand (Frequency)
	heatmap7 = sns.clustermap(FreqMinusID,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels7 = heatmap7.ax_heatmap.get_yticklabels()
	plt.setp(heatmap7.ax_heatmap.set_yticklabels(ylabels7,rotation=0))
	plt.setp(heatmap7.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap7.ax_heatmap.set_ylabel('{0} Elements'.format(len(FreqMinusID.index)),size=10))
	plt.setp(heatmap7.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap7.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap7.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap7.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap7.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap7.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap7.ax_heatmap.set_title('Methylation Frequency on Minus Strand for Elements',size=12))

	sns.despine()
	pp.savefig()
	print 'Plotted methylation frequency for element x position , element'
	
	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap8 = sns.clustermap(ranFreqPlusID,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels8 = heatmap8.ax_heatmap.get_yticklabels()
	plt.setp(heatmap8.ax_heatmap.set_yticks([]))
	plt.setp(heatmap8.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap8.ax_heatmap.set_ylabel('{0} Elements'.format(len(ranFreqPlusID.index)),size=10))
	plt.setp(heatmap8.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap8.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap8.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap8.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap8.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap8.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap8.ax_heatmap.set_title('Methylation Frequency on Plus Strand for Random Regions',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on neg strand (Frequency)
	heatmap9 = sns.clustermap(ranFreqMinusID,cmap='RdPu',xticklabels=50,col_cluster=False)
	ylabels9 = heatmap9.ax_heatmap.get_yticklabels()
	plt.setp(heatmap9.ax_heatmap.set_yticks([]))
	plt.setp(heatmap9.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap9.ax_heatmap.set_ylabel('{0} Elements'.format(len(ranFreqMinusID.index)),size=10))
	plt.setp(heatmap9.ax_heatmap.set_xlabel('Position',size=10))
	plt.setp(heatmap9.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap9.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap9.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap9.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5))
	plt.setp(heatmap9.ax_heatmap.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5))
	plt.setp(heatmap9.ax_heatmap.set_title('Methylation Frequency on Minus StrandStrand for Random Regions',size=12))

	sns.despine()
	pp.savefig()
	print 'Plotted methylation frequency for element x position , random regions'

	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap10 = sns.clustermap(XPlus,cmap='RdPu')
	ylabels10 = heatmap10.ax_heatmap.get_yticklabels()
	plt.setp(heatmap10.ax_heatmap.set_yticklabels(ylabels10,rotation=0))
	plt.setp(heatmap10.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap10.ax_heatmap.set_ylabel('{0} Elements'.format(len(FreqPlusID.index)),size=10))
	plt.setp(heatmap10.ax_heatmap.set_xlabel('Sample',size=10))
	plt.setp(heatmap10.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap10.ax_heatmap.set_title('Methylation Frequency on Plus Strand for Elements',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on neg strand (Frequency)
	heatmap11 = sns.clustermap(XMinus,cmap='RdPu')
	ylabels11 = heatmap11.ax_heatmap.get_yticklabels()
	plt.setp(heatmap11.ax_heatmap.set_yticklabels(ylabels11,rotation=0))
	plt.setp(heatmap11.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap11.ax_heatmap.set_ylabel('{0} Elements'.format(len(FreqMinusID.index)),size=10))
	plt.setp(heatmap11.ax_heatmap.set_xlabel('Sample',size=10))
	plt.setp(heatmap11.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap11.ax_heatmap.set_title('Methylation Frequency on Minus Strand for Elements',size=12))
	
	sns.despine()
	pp.savefig()
	print 'Plotted methylation frequency for element x tissue type , element'

	# Make heatmap for # methylation on pos strand (Frequency)
	heatmap12 = sns.clustermap(ranXPlus,cmap='RdPu')
	ylabels12 = heatmap12.ax_heatmap.get_yticklabels()
	plt.setp(heatmap12.ax_heatmap.set_yticks([]))
	plt.setp(heatmap12.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap12.ax_heatmap.set_ylabel('{0} Elements'.format(len(ranFreqPlusID.index)),size=10))
	plt.setp(heatmap12.ax_heatmap.set_xlabel('Sample',size=10))
	plt.setp(heatmap12.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap12.ax_heatmap.set_title('Methylation Frequency on Plus Strand for Random Regions',size=12))

	sns.despine()
	pp.savefig()
	
	# Make heatmap for # methylation on neg strand (Frequency)
	heatmap13 = sns.clustermap(ranXMinus,cmap='RdPu')
	ylabels13 = heatmap13.ax_heatmap.get_yticklabels()
	plt.setp(heatmap13.ax_heatmap.set_yticks([]))
	plt.setp(heatmap13.ax_heatmap.yaxis.tick_right())
	plt.setp(heatmap13.ax_heatmap.set_ylabel('{0} Elements'.format(len(ranFreqMinusID.index)),size=10))
	plt.setp(heatmap13.ax_heatmap.set_xlabel('Sample',size=10))
	plt.setp(heatmap13.ax_heatmap.tick_params(labelsize=10))
	plt.setp(heatmap13.ax_heatmap.set_title('Methylation Frequency on Minus Strand for Random Regions',size=12))
	
	sns.despine()
	pp.savefig()
	print 'Plotted methylation frequency for element x position , random regions'

# 	#put the index in a list
# 	UCEindex = ATelement.T.index.tolist()
# 	RANindex = ranATelement.T.index.tolist()
# 	
# 	reorder index based on clustering
# 	ATsorted = [UCEindex[i] for i in ATOrdered]
# 	RANsorted = [RANindex[i] for i in ranATOrdered]
# 
# 	GraphTableLibrary.main(ATOrdered,ranATOrdered,'Cluster_{0}'.format(fileName))
# 	print 'Created table for re-ordered mean AT cluster data'

	sns.despine()
	pp.savefig()
	pp.close()

def main(dfWindow,ranWindow,pdMeth,rnMeth,names,fileName):
	print 'Running GraphClusterLibrary'
	graphCluster(dfWindow,ranWindow,pdMeth,rnMeth,names,fileName)

if __name__ == "__main__":
	main()