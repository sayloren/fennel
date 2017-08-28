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
from GraphFangLibrary import collect_sum_two_nucleotides
import GraphTableLibrary
import GlobalVariables

# Reduce each id/tissue x location to get a corresponding tissue/id list of associations
def collect_index_by_association(dataframe,yItem,xItem):
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
def collect_tissue_by_id_dataframe(dataframe,yItem,zItem):

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
def collect_methylation_by_index(dataframe,yItem):

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
def make_dictionary_for_colors(ATelement,huslPalette):
	ATQcutPosition = pd.qcut(ATelement.std(axis=1),q=8,labels=False)
	ATQcutElement = pd.qcut(ATelement.std(),q=8,labels=False)
	lutElement = dict(zip(ATQcutElement.unique(), huslPalette))
	elementColors = ATQcutElement.map(lutElement)
	lutPosition = dict(zip(ATQcutPosition.unique(), huslPalette))
	positionColors = ATQcutPosition.map(lutPosition)
	print 'Made dictionary for standard deviation'
	return elementColors,positionColors

# Make some graphs for fangs
def graph_cluster(dfWindow,ranWindow,pdMeth,rnMeth,names,fileName):

	plt.figure(figsize=(7,7))

	# Get group, mean and standard deviation for AT
	ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(ranWindow,names,'A','T')
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
	elementColors,positionColors = make_dictionary_for_colors(ATelement,huslPalette)
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
	ranelementColors,ranpositionColors = make_dictionary_for_colors(ranATelement,huslPalette)
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
	FreqPlusID,FreqMinusID = collect_methylation_by_index(pdMeth,'id')
	FreqPlusTis,FreqMinusTis = collect_methylation_by_index(pdMeth,'tissue')
	XPlus,XMinus = collect_tissue_by_id_dataframe(pdMeth,'id','tissue')

	ranFreqPlusID,ranFreqMinusID = collect_methylation_by_index(rnMeth,'id')
	ranFreqPlusTis,ranFreqMinusTis = collect_methylation_by_index(rnMeth,'tissue')
	ranXPlus,ranXMinus = collect_tissue_by_id_dataframe(rnMeth,'id','tissue')

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
	print 'Running graph_clusterLibrary'
	graph_cluster(dfWindow,ranWindow,pdMeth,rnMeth,names,fileName)

if __name__ == "__main__":
	main()

"""
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "{}"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright {yyyy} {name of copyright owner}

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
   """