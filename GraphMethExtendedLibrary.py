"""
Script to print out all the extended methylation graphs

Wren Saylor
August 31 2017

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
from GraphMethLibrary import methIndex

# Make Methylation graphs
def graph_methylation_extended(pdMeth,rnMeth,fileName):
	sns.set_style('ticks')
# 	info = str(fileName) + ', '+ str(len(pdMeth['id'].unique())) + ' - ' "UCES"
# 	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Methylation_Extended_{0}.pdf'.format(fileName))
	plt.figure(figsize=(7,7))

	# Various combinations to plot on heatmaps
	FreqPlusTis,FreqMinusTis = methIndex(pdMeth,'tissue','methFreq')
	FreqPlusID,FreqMinusID = methIndex(pdMeth,'id','methFreq')
	ranFreqPlusTis,ranFreqMinusTis = methIndex(rnMeth,'tissue','methFreq')
	ranFreqPlusID,ranFreqMinusID = methIndex(rnMeth,'id','methFreq')

	PerPlusTis,PerMinusTis = methIndex(pdMeth,'tissue','methPer')
	CovPlusTis,CovMinusTis = methIndex(pdMeth,'tissue','methCov')
	ranPerPlusTis,ranPerMinusTis = methIndex(rnMeth,'tissue','methPer')
	ranCovPlusTis,ranCovMinusTis = methIndex(rnMeth,'tissue','methCov')

	PerPlusID,PerMinusID = methIndex(pdMeth,'id','methPer')
	CovPlusID,CovMinusID = methIndex(pdMeth,'id','methCov')
	ranPerPlusID,ranPerMinusID = methIndex(rnMeth,'id','methPer')
	ranCovPlusID,ranCovMinusID = methIndex(rnMeth,'id','methCov')

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
# 	Make heatmap for # methylation on pos strand (Frequency)
	ax0 = plt.subplot(gs[0,:])
	heatmap0 = sns.heatmap(FreqPlusTis,cmap='RdPu',ax=ax0,xticklabels=100)
	ax0.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax0.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax0.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax0.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax0.set_ylabel('Sample',size=8)
	ax0.set_xlabel('Position',size=6)
	ax0.tick_params(labelsize=8)
	ylabels0 = FreqPlusTis.index
	ax0.set_yticklabels(ylabels0,minor=False,rotation=0)
	ax0.set_yticks(np.arange(FreqPlusTis.shape[0]) + 0.5, minor=False)
	ax0.set_title('Methylation Frequency on Plus Strand for Elements',size=8)
# 	Make heatmap for # methylation on pos strand (Frequency)
	ax1 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap1 = sns.heatmap(FreqMinusTis,cmap='RdPu',ax=ax1,xticklabels=100)
	ax1.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax1.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax1.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax1.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax1.set_ylabel('Sample',size=8)
	ax1.set_xlabel('Position',size=6)
	ax1.tick_params(labelsize=8)
	ylabels1 = FreqMinusTis.index
	ax1.set_yticklabels(ylabels1,minor=False,rotation=0)
	ax1.set_yticks(np.arange(FreqMinusTis.shape[0]) + 0.5, minor=False)
	ax1.set_title('Methylation Frequency on Minus Strand for Elements',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
# 	Make heatmap for # methylation on pos strand (Frequency)
	ax2 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap2 = sns.heatmap(FreqPlusID,cmap='RdPu',ax=ax2,xticklabels=100)
	ax2.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax2.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax2.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax2.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax2.set_ylabel('{0} Elements'.format(len(FreqPlusID.index)),size=8)
	ax2.set_xlabel('Position',size=6)
	ax2.tick_params(labelsize=8)
	ax2.set_yticklabels([])
	ax2.set_yticks(np.arange(FreqPlusID.shape[0]) + 0.5, minor=False)
	ax2.set_title('Methylation Frequency on Plus Strand for Elements',size=8)
# 	Make heatmap for # methylation on neg strand (Frequency)
	ax3 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap3 = sns.heatmap(FreqMinusID,cmap='RdPu',ax=ax3,xticklabels=100)
	ax3.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax3.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax3.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax3.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax3.set_ylabel('{0} Elements'.format(len(FreqMinusID.index)),size=8)
	ax3.set_xlabel('Position',size=6)
	ax3.tick_params(labelsize=8)
	ax3.set_yticklabels([])
	ax3.set_yticks(np.arange(FreqMinusID.shape[0]) + 0.5, minor=False)
	ax3.set_title('Methylation Frequency on Minus Strand for Elements',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
# 	Make heatmap for # methylation on pos strand (Frequency)
	ax4 = plt.subplot(gs[0,:])
	heatmap4 = sns.heatmap(ranFreqPlusTis,cmap='RdPu',ax=ax4,xticklabels=100)#cbar_ax=cbar5_ax,vmin=0,vmax=5
	ax4.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax4.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax4.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax4.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax4.set_ylabel('Sample',size=8)
	ax4.set_xlabel('Position',size=6)
	ax4.tick_params(labelsize=8)
	ylabels4 = ranFreqPlusTis.index
	ax4.set_yticklabels(ylabels4,minor=False,rotation=0)
	ax4.set_yticks(np.arange(ranFreqPlusTis.shape[0]) + 0.5, minor=False)
	ax4.set_title('Methylation Frequency on Plus Strand for Random Regions',size=8)
# 	Make heatmap for # methylation on pos strand (Frequency)
	ax5 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap5 = sns.heatmap(ranFreqMinusTis,cmap='RdPu',ax=ax5,xticklabels=100)#cbar_ax=cbar5_ax,vmin=0,vmax=5
	ax5.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax5.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax5.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax5.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax5.set_ylabel('Sample',size=8)
	ax5.set_xlabel('Position',size=6)
	ax5.tick_params(labelsize=8)
	ylabels5 = ranFreqMinusTis.index
	ax5.set_yticklabels(ylabels5,minor=False,rotation=0)
	ax5.set_yticks(np.arange(ranFreqMinusTis.shape[0]) + 0.5, minor=False)
	ax5.set_title('Methylation Frequency on Minus Strand for Random Regions',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
# 	Make heatmap for # methylation on pos strand (Frequency)
	ax6 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap6 = sns.heatmap(ranFreqPlusID,cmap='RdPu',ax=ax6,xticklabels=100)#,vmin=0,vmax=5
	ax6.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax6.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax6.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax6.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax6.set_ylabel('{0} Elements'.format(len(ranFreqPlusID.index)),size=8)
	ax6.set_xlabel('Position',size=6)
	ax6.tick_params(labelsize=8)
	ax6.set_yticklabels([])
	ax6.set_yticks(np.arange(ranFreqPlusID.shape[0]) + 0.5, minor=False)
	ax6.set_title('Methylation Frequency on Plus Strand for Random Regions',size=8)
# 	Make heatmap for # methylation on neg strand (Frequency)
	ax7 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap7 = sns.heatmap(ranFreqMinusID,cmap='RdPu',ax=ax7,xticklabels=100)#,vmin=0,vmax=5
	ax7.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax7.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax7.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax7.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax7.set_ylabel('{0} Elements'.format(len(ranFreqMinusID.index)),size=8)
	ax7.set_xlabel('Position',size=6)
	ax7.tick_params(labelsize=8)
	ax7.set_yticklabels([])
	ax7.set_yticks(np.arange(ranFreqMinusID.shape[0]) + 0.5, minor=False)
	ax7.set_title('Methylation Frequency on Minus Strand for Random Regions',size=8)
	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	# Make heatmap for % methylation on pos strand
	# Should make the bottom threshold a different color
	ax8 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap8 = sns.heatmap(PerPlusTis,cmap='RdPu',ax=ax8,vmin=0,vmax=100,xticklabels=100)
	ax8.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax8.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax8.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax8.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax8.set_ylabel('Sample',size=8)
	ax8.set_xlabel('Position',size=6)
	ax8.tick_params(labelsize=8)
	ylabels8 = PerPlusTis.index
	ax8.set_yticklabels(ylabels8,minor=False,rotation=0)
	ax8.set_yticks(np.arange(PerPlusTis.shape[0]) + 0.5, minor=False)
	ax8.set_title('Methylation Percentage on Plus Strand for Element',size=8)
	# Make heatmap for % methylation on neg strand
	# Should make the bottom threshold a different color
	ax9 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap9 = sns.heatmap(PerMinusTis,cmap='RdPu',ax=ax9,vmin=0,vmax=100,xticklabels=100)
	ax9.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax9.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax9.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax9.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax9.set_ylabel('Sample',size=8)
	ax9.set_xlabel('Position',size=6)
	ax9.tick_params(labelsize=8)
	ylabels9 = PerMinusTis.index
	ax9.set_yticklabels(ylabels9,minor=False,rotation=0)
	ax9.set_yticks(np.arange(PerMinusTis.shape[0]) + 0.5, minor=False)
	ax9.set_title('Methylation Percentage on Minus Strand for Element',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	# Make heatmap for methylation coverage on pos strand, capped at 1000
	# Should make the bottom threshold a different color
	ax10 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap10 = sns.heatmap(CovPlusTis,cmap='RdPu',ax=ax10,vmin=0,vmax=100,xticklabels=100)
	ax10.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax10.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax10.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax10.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax10.set_ylabel('Sample',size=8)
	ax10.set_xlabel('Position',size=6)
	ax10.tick_params(labelsize=8)
	ylabels10 = CovPlusTis.index
	ax10.set_yticklabels(ylabels10,minor=False,rotation=0)
	ax10.set_yticks(np.arange(CovPlusTis.shape[0]) + 0.5, minor=False)
	ax10.set_title('Methylation Coverage for Plus Strand for Element',size=8)
	# Make heatmap for methylation coverage on neg strand, capped at 1000
	# Should make the bottom threshold a different color
	ax11 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap11 = sns.heatmap(CovMinusTis,cmap='RdPu',ax=ax11,vmin=0,vmax=100,xticklabels=100)
	ax11.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax11.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax11.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax11.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax11.set_ylabel('Sample',size=8)
	ax11.set_xlabel('Position',size=6)
	ax11.tick_params(labelsize=8)
	ylabels11 = CovMinusTis.index
	ax11.set_yticklabels(ylabels11,minor=False,rotation=0)
	ax11.set_yticks(np.arange(CovMinusTis.shape[0]) + 0.5, minor=False)
	ax11.set_title('Methylation Coverage for Minus Strand for Element',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	# Make heatmap for % methylation on pos strand
	# Should make the bottom threshold a different color
	ax12 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap12 = sns.heatmap(ranPerPlusTis,cmap='RdPu',ax=ax12,vmin=0,vmax=100,xticklabels=100)
	ax12.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax12.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax12.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax12.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax12.set_ylabel('Sample',size=8)
	ax12.set_xlabel('Position',size=6)
	ax12.tick_params(labelsize=8)
	ylabels12 = ranPerPlusTis.index
	ax12.set_yticklabels(ylabels12,minor=False,rotation=0)
	ax12.set_yticks(np.arange(ranPerPlusTis.shape[0]) + 0.5, minor=False)
	ax12.set_title('Methylation Percentage on Plus Strand for Random Regions',size=8)
	# Make heatmap for % methylation on neg strand
	# Should make the bottom threshold a different color
	ax13 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap13 = sns.heatmap(ranPerMinusTis,cmap='RdPu',ax=ax13,vmin=0,vmax=100,xticklabels=100)
	ax13.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax13.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax13.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax13.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax13.set_ylabel('Sample',size=8)
	ax13.set_xlabel('Position',size=6)
	ax13.tick_params(labelsize=8)
	ylabels13 = ranPerMinusTis.index
	ax13.set_yticklabels(ylabels13,minor=False,rotation=0)
	ax13.set_yticks(np.arange(ranPerMinusTis.shape[0]) + 0.5, minor=False)
	ax13.set_title('Methylation Percentage on Minus Strand for Random Regions',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	# Make heatmap for methylation coverage on pos strand, capped at 1000
	# Should make the bottom threshold a different color
	ax14 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap14 = sns.heatmap(ranCovPlusTis,cmap='RdPu',ax=ax14,vmin=0,vmax=100,xticklabels=100)
	ax14.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax14.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax14.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax14.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax14.set_ylabel('Sample',size=8)
	ax14.set_xlabel('Position',size=6)
	ax14.tick_params(labelsize=8)
	ylabels14 = ranCovPlusTis.index
	ax14.set_yticklabels(ylabels14,minor=False,rotation=0)
	ax14.set_yticks(np.arange(ranCovPlusTis.shape[0]) + 0.5, minor=False)
	ax14.set_title('Methylation Coverage for Plus Strand for Random Regions',size=8)
	# Make heatmap for methylation coverage on neg strand, capped at 1000
	# Should make the bottom threshold a different color
	ax15 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap15 = sns.heatmap(ranCovMinusTis,cmap='RdPu',ax=ax15,vmin=0,vmax=100,xticklabels=100)
	ax15.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax15.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax15.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax15.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax15.set_ylabel('Sample',size=8)
	ax15.set_xlabel('Position',size=6)
	ax15.tick_params(labelsize=8)
	ylabels15 = ranCovMinusTis.index
	ax15.set_yticklabels(ylabels15,minor=False,rotation=0)
	ax15.set_yticks(np.arange(ranCovMinusTis.shape[0]) + 0.5, minor=False)
	ax15.set_title('Methylation Coverage for Minus Strand for Random Regions',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	# Make heatmap for % methylation on pos strand
	# Should make the bottom threshold a different color
	ax16 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap16 = sns.heatmap(PerPlusID,cmap='RdPu',ax=ax16,vmin=0,vmax=100,xticklabels=100)
	ax16.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax16.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax16.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax16.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax16.set_ylabel('{0} Elements'.format(len(PerPlusID.index)),size=8)
	ax16.set_xlabel('Position',size=6)
	ax16.tick_params(labelsize=8)
	ax16.set_yticklabels([])
	ax16.set_title('Methylatoin Percentage on Plus Strand for Element',size=8)
	# Make heatmap for % methylation on neg strand
	# Should make the bottom threshold a different color
	ax17 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap17 = sns.heatmap(PerMinusID,cmap='RdPu',ax=ax17,vmin=0,vmax=100,xticklabels=100)
	ax17.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax17.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax17.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax17.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax17.set_ylabel('{0} Elements'.format(len(PerMinusID.index)),size=8)
	ax17.set_xlabel('Position',size=6)
	ax17.tick_params(labelsize=8)
	ax17.set_yticklabels([])
	ax17.set_title('Methylatoin Percentage on Minus Strand for Element',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	# Make heatmap for methylation coverage on pos strand, capped at 1000
	# Should make the bottom threshold a different color
	ax18 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap18 = sns.heatmap(CovPlusID,cmap='RdPu',ax=ax18,xticklabels=100,vmin=0,vmax=100)
	ax18.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax18.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax18.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax18.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax18.set_ylabel('{0} Elements'.format(len(CovPlusID.index)),size=8)
	ax18.set_xlabel('Position',size=6)
	ax18.tick_params(labelsize=8)
	ax18.set_yticklabels([])
	ax18.set_title('Methylation Coverage on Plus Strand for Element',size=8)
	# Make heatmap for methylation coverage on neg strand, capped at 1000
	# Should make the bottom threshold a different color
	ax19 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap7 = sns.heatmap(CovMinusID,cmap='RdPu',ax=ax19,xticklabels=100,vmin=0,vmax=100)
	ax19.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax19.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax19.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax19.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax19.set_ylabel('{0} Elements'.format(len(CovMinusID.index)),size=8)
	ax19.set_xlabel('Position',size=6)
	ax19.tick_params(labelsize=8)
	ax19.set_yticklabels([])
	ax19.set_title('Methylation Coverage on Minus Strand for Element',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	# Make heatmap for % methylation on pos strand
	# Should make the bottom threshold a different color
	ax20 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap20 = sns.heatmap(ranPerPlusID,cmap='RdPu',ax=ax20,vmin=0,vmax=100,xticklabels=100)
	ax20.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax20.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax20.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax20.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax20.set_ylabel('{0} Elements'.format(len(ranPerPlusID.index)),size=8)
	ax20.set_xlabel('Position',size=6)
	ax20.tick_params(labelsize=8)
	ax20.set_yticklabels([])
	ax20.set_title('Methylatoin Percentage on Plus Strand for Random Regions',size=8)
	# Make heatmap for % methylation on neg strand
	# Should make the bottom threshold a different color
	ax21 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap11 = sns.heatmap(ranPerMinusID,cmap='RdPu',ax=ax21,vmin=0,vmax=100,xticklabels=100)
	ax21.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax21.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax21.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax21.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax21.set_ylabel('{0} Elements'.format(len(ranPerMinusID.index)),size=8)
	ax21.set_xlabel('Position',size=6)
	ax21.tick_params(labelsize=8)
	ax21.set_yticklabels([])
	ax21.set_title('Methylatoin Percentage on Minus Strand for Random Regions',size=8)
	sns.despine()
	pp.savefig()

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	# Make heatmap for methylation coverage on pos strand, capped at 1000
	# Should make the bottom threshold a different color
	ax22 = plt.subplot(gs[0,:],sharex=ax0)
	heatmap22 = sns.heatmap(ranCovPlusID,cmap='RdPu',ax=ax22,xticklabels=100,vmin=0,vmax=100)
	ax22.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax22.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax22.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax22.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax22.set_ylabel('{0} Elements'.format(len(ranCovPlusID.index)),size=8)
	ax22.set_xlabel('Position',size=6)
	ax22.tick_params(labelsize=8)
	ax22.set_yticklabels([])
	ax22.set_title('Methylation Coverage on Plus Strand for Random Regions',size=8)
	# Make heatmap for methylation coverage on neg strand, capped at 1000
	# Should make the bottom threshold a different color
	ax23 = plt.subplot(gs[1,:],sharex=ax0)
	heatmap23 = sns.heatmap(ranCovMinusID,cmap='RdPu',ax=ax23,xticklabels=100,vmin=0,vmax=100)
	ax23.axvline(x=GlobalVariables.plotLineLocationOneFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax23.axvline(x=GlobalVariables.plotLineLocationTwoFull,linewidth=.05,linestyle='dashed',color='#5fc85b',alpha=0.5)
	ax23.axvline(x=GlobalVariables.plotLineLocationThreeFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax23.axvline(x=GlobalVariables.plotLineLocationFourFull,linewidth=.05,linestyle='dashed',color='#96c85b',alpha=0.5)
	ax23.set_ylabel('{0} Elements'.format(len(ranCovMinusID.index)),size=8)
	ax23.set_xlabel('Position',size=6)
	ax23.tick_params(labelsize=8)
	ax23.set_yticklabels([])
	ax23.set_title('Methylation Coverage on Minus Strand for Random Regions',size=8)


	sns.despine()
	pp.savefig()
	pp.close()
	print 'Plotted methylation frequency across elements'

def main(pdMeth,rnMeth,fileName):
	graph_methylation_extended(pdMeth,rnMeth,fileName)

if __name__ == "__main__":
	main()
