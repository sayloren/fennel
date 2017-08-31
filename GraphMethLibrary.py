"""
Script to graph the methylation analyses

Wren Saylor
July 5 2017

Copyright 2017 Harvard University, Wu Lab

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
import GlobalVariables

# Transform the Frequency, Percentage and Coverage data into graphable data frames
def methIndex(dataframe,yItem,zItem): 
	
	# x item is methLoc, y item is either tissue or id, z item is coverage, percentage, or frequency
	new_index = range(0,GlobalVariables.num)
	
	# Separate by strand
	PlusMeth = dataframe.loc[dataframe['Cytosine'] == 'C']
	MinusMeth = dataframe.loc[dataframe['Cytosine'] == 'G']
	
	# Subset just columns to use
	PlussubMeth = PlusMeth[['methLoc',yItem,zItem]]
	MinussubMeth = MinusMeth[['methLoc',yItem,zItem]]
	
	# Sort ascending, in order to only use the highest value with keep = last
	PlussortMeth = PlussubMeth.sort_values(['methLoc'],ascending=True).drop_duplicates(['methLoc',yItem,zItem],keep='last')
	MinussortMeth = MinussubMeth.sort_values(['methLoc'],ascending=True).drop_duplicates(['methLoc',yItem,zItem],keep='last')

	# Pivot the data frame so that each tissue/cell type is a column
	PluspivotMeth = pd.pivot_table(PlussortMeth,index='methLoc',columns=yItem,values=zItem)
	MinuspivotMeth = pd.pivot_table(MinussortMeth,index='methLoc',columns=yItem,values=zItem)
	
	PluspivotMeth.columns.name = None
	MinuspivotMeth.columns.name = None
	
	# Give new index, using the methLocations
	PlusindexMeth = PluspivotMeth.reindex(new_index,fill_value=0)
	MinusindexMeth = MinuspivotMeth.reindex(new_index,fill_value=0)

	# Fill in missing index with 0s
	PlusindexMeth.fillna('0',inplace=True)
	MinusindexMeth.fillna('0',inplace=True)
	
	# Remove the index column name
	PlusindexMeth.index.name = None
	MinusindexMeth.index.name = None
	
	# Transpose the data frame for easy input into the heatamp
	PlustransMeth = PlusindexMeth.T
	MinustransMeth = MinusindexMeth.T
	
	PlustransMeth = PlustransMeth[PlustransMeth.columns].astype(float)
	MinustransMeth = MinustransMeth[MinustransMeth.columns].astype(float)
	
	print 'Converted {0} by {1} into data frame'.format(yItem,zItem)
	
	return PlustransMeth, MinustransMeth

# Make Methylation graphs
def graphMeth(pdMeth,rnMeth,fileName):
	sns.set_style('ticks')
	info = str(fileName) + ', '+ str(len(pdMeth['id'].unique())) + ' - ' "UCES"
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Methylation_{0}.pdf'.format(fileName))
	plt.figure(figsize=(7,7))

	# Various combinations to plot on heatmaps
	FreqPlusTis,FreqMinusTis = methIndex(pdMeth,'tissue','methFreq')
	FreqPlusID,FreqMinusID = methIndex(pdMeth,'id','methFreq')
	ranFreqPlusTis,ranFreqMinusTis = methIndex(rnMeth,'tissue','methFreq')
	ranFreqPlusID,ranFreqMinusID = methIndex(rnMeth,'id','methFreq')

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	nostrandFreqTis = FreqPlusTis.add(FreqMinusTis,fill_value=0)
	nostrandranFreqTis = ranFreqPlusTis.add(ranFreqMinusTis,fill_value=0)
	nostrandFreqID = FreqPlusID.add(FreqMinusID,fill_value=0)
	nostrandranFreqID = ranFreqPlusID.add(ranFreqMinusID,fill_value=0)

# 	Make heatmap for # methylation on pos strand (Frequency)
	ax0 = plt.subplot(gs[0,:])
	heatmap0 = sns.heatmap(nostrandFreqTis,cmap='RdPu',ax=ax0,xticklabels=100)
	ax0.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax0.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax0.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax0.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax0.set_ylabel('Sample',size=8)
	ax0.set_xlabel('Position',size=6)
	ax0.tick_params(labelsize=8)
	ylabels0 = nostrandFreqTis.index
	ax0.set_yticklabels(ylabels0,minor=False,rotation=0)
	ax0.set_yticks(np.arange(nostrandFreqTis.shape[0]) + 0.5, minor=False)
	ax0.set_title('Methylation Frequency for Elements',size=8)

# 	Make heatmap for # methylation on pos strand (Frequency)
	ax1 = plt.subplot(gs[1,:])
	heatmap1 = sns.heatmap(nostrandranFreqTis,cmap='RdPu',ax=ax1,xticklabels=100)
	ax1.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax1.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax1.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax1.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax1.set_ylabel('Sample',size=8)
	ax1.set_xlabel('Position',size=6)
	ax1.tick_params(labelsize=8)
	ylabels1 = nostrandranFreqTis.index
	ax1.set_yticklabels(ylabels1,minor=False,rotation=0)
	ax1.set_yticks(np.arange(nostrandranFreqTis.shape[0]) + 0.5, minor=False)
	ax1.set_title('Methylation Frequency for Random Regions',size=8)

	sns.despine()
	pp.savefig()
	print 'Plotted methylation frequency across tissue types'

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

# 	Make heatmap for # methylation on pos strand (Frequency)
	ax2 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap2 = sns.heatmap(nostrandFreqID,cmap='RdPu',ax=ax2,xticklabels=100)
	ax2.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax2.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax2.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax2.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax2.set_ylabel('{0} Elements'.format(len(nostrandFreqID.index)),size=8)
	ax2.set_xlabel('Position',size=6)
	ax2.tick_params(labelsize=8)
	ax2.set_yticklabels([])
	ax2.set_yticks(np.arange(nostrandFreqID.shape[0]) + 0.5, minor=False)
	ax2.set_title('Methylation Frequency for Elements',size=8)

# 	Make heatmap for # methylation on pos strand (Frequency)
	ax3 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap3 = sns.heatmap(nostrandranFreqID,cmap='RdPu',ax=ax3,xticklabels=100)
	ax3.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax3.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax3.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax3.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax3.set_ylabel('{0} Elements'.format(len(nostrandranFreqID.index)),size=8)
	ax3.set_xlabel('Position',size=6)
	ax3.tick_params(labelsize=8)
	ax3.set_yticklabels([])
	ax3.set_yticks(np.arange(nostrandranFreqID.shape[0]) + 0.5, minor=False)
	ax3.set_title('Methylation Frequency for Random Regions',size=8)

	sns.despine()
	pp.savefig()
	pp.close()
	print 'Plotted methylation frequency across elements'

def main(pdMeth,rnMeth,fileName):
	print 'Running GraphMethylationLibrary'
	graphMeth(pdMeth,rnMeth,fileName)

if __name__ == "__main__":
	main()
