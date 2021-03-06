"""
Script to graph the Fang analyses

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
from scipy.interpolate import splrep, splev
import scipy.stats as ss
from scipy.stats import mstats
import seaborn as sns
import GlobalVariables

# Get group, mean and standard deviation for AT
def collect_sum_two_nucleotides(dfWindow,names,nuc1,nuc2):
	ATNames = [names.index(i) for i in names if nuc1 in i or nuc2 in i]
	ATDataFrames = [dfWindow[i] for i in ATNames]
	ATconcat = pd.concat(ATDataFrames,axis=1)
	ATgroup = ATconcat.groupby(ATconcat.columns,axis=1).sum()
	ATmean = ATgroup.mean()
	ATstd = ATgroup.std()
	return ATgroup, ATmean, ATstd

# Make some graphs for fangs
def graph_element_line_means(dfWindow,names,ranWindow,fileName):

	# Get group, mean and standard deviation for AT
	ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(ranWindow,names,'A','T')
	
	# Title info
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	gs = gridspec.GridSpec(3,1,height_ratios=[3,1,1])
	gs.update(hspace=.8) # setting the space between the graphs
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Fangs_{0}.pdf'.format(fileName))
	plt.figure(figsize=(7,7))

	sns.set_palette("husl",n_colors=8)#(len(nucLine)*2)

	# Plot the mean AT content with a std of 1
# 	StartMean = ATgroup.loc[:,(GlobalVariables.plotLineLocationThree-GlobalVariables.window):GlobalVariables.plotLineLocationOne].mean() # the boundary and 50inward
# 	StopMean = ATgroup.loc[:,(GlobalVariables.plotLineLocationTwo):(GlobalVariables.plotLineLocationFour-GlobalVariables.window)].mean() # the boundary and 50inward
# 	wilcoxPSRMean = ss.wilcoxon(StartMean,StopMean)
	ax0 = plt.subplot(gs[0])
	ax0.plot(GlobalVariables.fillX,ATmean,linewidth=1,label='AT element')
	ax0.plot(GlobalVariables.fillX,ranATmean,linewidth=1,label='AT random')
	ax0.fill_between(GlobalVariables.fillX,ATmean+ATstd,ATmean-ATstd,label='',alpha=0.2)
	ax0.fill_between(GlobalVariables.fillX,ranATmean+ranATstd,ranATmean-ranATstd,label='',alpha=0.2)
	ax0.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.hlines(y=86,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
	ax0.text(32,85,'11bp sliding window',size=6)
# 	ax0.text(20,90,'Wilcox Signed Rank P-value {:0.1e}'.format(wilcoxPSRMean[1]),size=6,clip_on=False)
	ax0.set_ylabel('% AT Content',size=8)
	ax0.set_xlabel('Position',size=6)
	ax0.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax0.set_title('Mean AT Content With Standard Deviation',size=8)
	ax0.set_yticks(ax0.get_yticks()[::2])
	plt.xlim(0,GlobalVariables.num)

	# Plot the std = 1
	ax1 = plt.subplot(gs[1],sharex=ax0)
	ax1.plot(GlobalVariables.fillX,ATstd,linewidth=1,label='AT element')
	ax1.plot(GlobalVariables.fillX,ranATstd,linewidth=1,label='AT random')
	ax1.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.set_yticks(ax1.get_yticks()[::2])
	ax1.set_xlabel('Position',size=6)
	ax1.set_ylabel('SD',size=8)
	ax1.set_title('Standard Deviation',size=8)
	plt.setp(ax1.get_xticklabels(), visible=True)
	ax1.legend(loc=0,fontsize=5,labelspacing=0.05)
	
	# Significances tests for SD populations
	ATuceRegion = ATgroup.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
	ATranuceRegion = ranATgroup.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
	ATbothStream = ATgroup.iloc[:,np.r_[0:GlobalVariables.plotLineLocationThree,GlobalVariables.plotLineLocationFour:(GlobalVariables.num-GlobalVariables.window)]].std()
	ATranbothStream = ranATgroup.iloc[:,np.r_[0:GlobalVariables.plotLineLocationThree,GlobalVariables.plotLineLocationFour:(GlobalVariables.num-GlobalVariables.window)]].std()
	# Kruskal-Wallis test
# 	kruskalSD = mstats.kruskalwallis(bothStream,uceRegion)
# 	ax2.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)
	ax2 = plt.subplot(gs[2])
	ax2.hist(ATuceRegion,35,linewidth=0.3,label='Element',alpha=0.4)
	ax2.hist(ATranuceRegion,35,linewidth=0.3,label='Random',alpha=0.4)
	ax2.hist(ATbothStream,35,linewidth=0.3,label='Surrounding Element Regions',alpha=0.4)
	ax2.hist(ATranbothStream,35,linewidth=0.3,label='Surrounding Random Regions',alpha=0.4)
	ax2.set_yticks(ax2.get_yticks()[::2])
	ax2.set_ylabel('Frequency',size=8)
	ax2.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax2.set_xlabel('Standard Deviation Value',size=8)

	sns.despine()
	plt.savefig(pp, format='pdf')
	print 'Plotted the mean AT content and standard deviation'

	# Separate out those with only a single nucleotide search
	SingleNames = [names.index(i) for i in names if len(i) == 1]
	SingleNamesVal = [names[i] for i in SingleNames]
	SingleDataFrames = [dfWindow[i] for i in SingleNames]
	ranSingleDataFrames = [ranWindow[i] for i in SingleNames]
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	
	# Nucleotides
	ax3 = plt.subplot(gs[0,:],sharex=ax0)
	for dfNuc,lNuc in zip(SingleDataFrames,SingleNamesVal):
		ax3.plot(GlobalVariables.fillX,dfNuc.mean(),linewidth=1,label=lNuc)
	ax3.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax3.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax3.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax3.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax3.set_yticks(ax3.get_yticks()[::2])
	ax3.set_ylabel('% Nucleotide Content',size=8)
	ax3.set_xlabel('Position',size=6)
	ax3.set_title('Mean Nucleotide Content for Element',size=8)
	ax3.legend(loc=0,fontsize=5,labelspacing=0.05)
	
	# Plot SD for Nucleotides
	ax4 = plt.subplot(gs[1,:])
	for dfNuc,lNuc in zip(SingleDataFrames,SingleNamesVal):
		elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
		ax4.hist(elRegion,35,linewidth=0.3,alpha=0.5,label='{0}-{1}'.format(lNuc,'element'))
	ax4.set_yticks(ax4.get_yticks()[::2])
	ax4.set_ylabel('Frequency',size=8)
	ax4.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax4.set_xlabel('Standard Deviation Value',size=8)
	ax4.set_title('Standard Deviation for Nucleotides within Element',size=8)
	
	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)
	
	# Nucleotides
	ax5 = plt.subplot(gs[0,:],sharex=ax0)
	for dfNuc,lNuc in zip(ranSingleDataFrames,SingleNamesVal):
		ax5.plot(GlobalVariables.fillX,dfNuc.mean(),linewidth=1,label=lNuc)
	ax5.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax5.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax5.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax5.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax5.set_yticks(ax5.get_yticks()[::2])
	ax5.set_ylabel('% Nucleotide Content',size=8)
	ax5.set_xlabel('Position',size=6)
	ax5.set_title('Mean Nucleotide Content for Random Region',size=8)
	ax5.legend(loc=0,fontsize=5,labelspacing=0.05)
	
	# Plot SD for Nucleotides
	ax6 = plt.subplot(gs[1,:])
	for dfNuc,lNuc in zip(ranSingleDataFrames,SingleNamesVal):
		elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
		ax6.hist(elRegion,35,linewidth=0.3,alpha=0.5,label='{0}-{1}'.format(lNuc,'element'))
	ax6.set_yticks(ax6.get_yticks()[::2])
	ax6.set_ylabel('Frequency',size=8)
	ax6.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax6.set_xlabel('Standard Deviation Value',size=8)
	ax6.set_title('Standard Deviation for Nucleotides within Random Region',size=8)
	
	print 'Plotted nucleotide means'

	# For any search stings that are dinucleotides
	if any(len(i) == 2 for i in names):
		
		sns.despine()
		pp.savefig()
		
		# Plot settings
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)
		
		# Separate out those with only a double nucleotide search
		DoubleNames = [names.index(i) for i in names if len(i) == 2]
		DoubleNamesVal = [names[i] for i in DoubleNames]
		DoubleDataFrames = [dfWindow[i] for i in DoubleNames]
		ranDoubleDataFrames = [ranWindow[i] for i in DoubleNames]
		
		# Might still want to return the actual CpN location for how many are methylated
		ax7 = plt.subplot(gs[0,:],sharex=ax0)
		for dfNuc,lNuc in zip(DoubleDataFrames,DoubleNamesVal):
			ax7.plot(GlobalVariables.fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax7.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax7.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax7.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax7.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax7.set_title('Mean Nucleotide String Content for Element',size=8)
		ax7.set_ylabel('% Nucleotide String Content',size=8)
		ax7.set_xlabel('Position',size=6)
		ax7.set_yticks(ax7.get_yticks()[::2])
		ax7.legend(loc=0,fontsize=5,labelspacing=0.05)

		ax8 = plt.subplot(gs[1,:])
		for dfNuc,lNuc in zip(DoubleDataFrames,DoubleNamesVal):
			elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
			ax8.hist(elRegion,25,linewidth=0.3,label=lNuc,alpha=0.5)
		ax8.set_yticks(ax8.get_yticks()[::2])
		ax8.set_ylabel('Frequency',size=8)
		ax8.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax8.set_xlabel('Standard Deviation Value',size=8)
		ax8.set_title('Standard Deviation for Nucleotide String within Element',size=8)

		sns.despine()
		pp.savefig()
		
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)

		ax9 = plt.subplot(gs[0,:],sharex=ax0)
		for dfNuc,lNuc in zip(ranDoubleDataFrames,DoubleNamesVal):
			ax9.plot(GlobalVariables.fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax9.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax9.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax9.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax9.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax9.set_title('Mean Nucleotide String Content for Random Region',size=8)
		ax9.set_ylabel('% Nucleotide String Content',size=8)
		ax9.set_xlabel('Position',size=6)
		ax9.set_yticks(ax9.get_yticks()[::2])
		ax9.legend(loc=0,fontsize=5,labelspacing=0.05)

		ax10 = plt.subplot(gs[1,:])
		for dfNuc,lNuc in zip(ranDoubleDataFrames,DoubleNamesVal):
			elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
			ax10.hist(elRegion,25,linewidth=0.3,label=lNuc,alpha=0.5)
		ax10.set_yticks(ax10.get_yticks()[::2])
		ax10.set_ylabel('Frequency',size=8)
		ax10.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax10.set_xlabel('Standard Deviation Value',size=8)
		ax10.set_title('Standard Deviation for Nucleotide String within Random Region',size=8)
		
		print 'Plotted dinucleotide means'

	# For any search strings 3 bp long
	if any(len(i) == 3 for i in names):
	
		sns.despine()
		pp.savefig()
		
		# Plot settings
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)
		
		# Separate out those with only a multiple nucleotide search
		TriNames = [names.index(i) for i in names if len(i) == 3]
		TriNamesVal = [names[i] for i in TriNames]
		TriDataFrames = [dfWindow[i] for i in TriNames]
		ranTriDataFrames = [ranWindow[i] for i in TriNames]

		# Plot the MultiNucleotide Sequences
		ax11 = plt.subplot(gs[0],sharex=ax0)
		for dfNuc,lNuc in zip(TriDataFrames,TriNamesVal):
			ax11.plot(GlobalVariables.fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax11.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax11.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax11.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax11.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax11.set_title('Mean Trinucleotide Content for Element',size=8)
		ax11.set_ylabel('% Trinucleotide Content',size=8)
		ax11.set_xlabel('Position',size=6)
		ax11.set_yticks(ax11.get_yticks()[::2])
		ax11.legend(loc=0,fontsize=5,labelspacing=0.05)

# 		Plot SD for MultiNucleotide Sequences
		ax12 = plt.subplot(gs[1])
		for dfNuc,lNuc in zip(TriDataFrames,TriNamesVal):
			elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
			ax12.hist(elRegion,20,linewidth=0.3,label=lNuc,alpha=0.5)
		ax12.set_yticks(ax12.get_yticks()[::2])
		ax12.set_ylabel('Frequency',size=8)
		ax12.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax12.set_xlabel('Standard Deviation Value',size=8)
		ax12.set_title('Standard Deviation for Trinucleotides within Element',size=8)
		
		sns.despine()
		pp.savefig()
		
		# Plot settings
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)

		# Plot the MultiNucleotide Sequences
		ax13 = plt.subplot(gs[0],sharex=ax0)
		for dfNuc,lNuc in zip(ranTriDataFrames,TriNamesVal):
			ax13.plot(GlobalVariables.fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax13.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax13.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax13.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax13.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax13.set_title('Mean Trinucleotide Content for Random Region',size=8)
		ax13.set_ylabel('% Trinucleotide Content',size=8)
		ax13.set_xlabel('Position',size=6)
		ax13.set_yticks(ax13.get_yticks()[::2])
		ax13.legend(loc=0,fontsize=5,labelspacing=0.05)

# 		Plot SD for MultiNucleotide Sequences
		ax14 = plt.subplot(gs[1])
		for dfNuc,lNuc in zip(ranTriDataFrames,TriNamesVal):
			elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
			ax14.hist(elRegion,20,linewidth=0.3,label=lNuc,alpha=0.5)
		ax14.set_yticks(ax14.get_yticks()[::2])
		ax14.set_ylabel('Frequency',size=8)
		ax14.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax14.set_xlabel('Standard Deviation Value',size=8)
		ax14.set_title('Standard Deviation for Trinucleotides within Random Region',size=8)

		print 'Plotted trinucleotide means'
		
	# For any search strings longer than 3
	if any(len(i) > 3 for i in names):
	
		sns.despine()
		pp.savefig()
		
		# Plot settings
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)
		
		# Separate out those with only a multiple nucleotide search
		MultiNames = [names.index(i) for i in names if len(i) > 3]
		MultiNamesVal = [names[i] for i in MultiNames]
		MultiDataFrames = [dfWindow[i] for i in MultiNames]
		ranMultiDataFrames = [ranWindow[i] for i in MultiNames]

		# Plot the MultiNucleotide Sequences
		ax15 = plt.subplot(gs[0],sharex=ax0)
		for dfNuc,lNuc in zip(MultiDataFrames,MultiNamesVal):
			ax15.plot(GlobalVariables.fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax15.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax15.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax15.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax15.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax15.set_title('Mean Polynucleotide Content for Element',size=8)
		ax15.set_ylabel('% Polynucleotide Content',size=8)
		ax15.set_xlabel('Position',size=6)
		ax15.set_yticks(ax15.get_yticks()[::2])
		ax15.legend(loc=0,fontsize=5,labelspacing=0.05)

# 		Plot SD for MultiNucleotide Sequences
		ax16 = plt.subplot(gs[1])
		for dfNuc,lNuc in zip(MultiDataFrames,MultiNamesVal):
			elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
			ax16.hist(elRegion,20,linewidth=0.3,label=lNuc,alpha=0.5)
		ax16.set_yticks(ax16.get_yticks()[::2])
		ax16.set_ylabel('Frequency',size=8)
		ax16.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax16.set_xlabel('Standard Deviation Value',size=8)
		ax16.set_title('Standard Deviation for Polynucleotides within Element',size=8)
		
		sns.despine()
		pp.savefig()
		
		# Plot settings
		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
		gs.update(hspace=.5)

		# Plot the MultiNucleotide Sequences
		ax17 = plt.subplot(gs[0],sharex=ax0)
		for dfNuc,lNuc in zip(ranMultiDataFrames,MultiNamesVal):
			ax17.plot(GlobalVariables.fillX,dfNuc.mean(),linewidth=1,label=lNuc)
		ax17.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax17.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
		ax17.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax17.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
		ax17.set_title('Mean Polynucleotide Content for Random Region',size=8)
		ax17.set_ylabel('% Polynucleotide Content',size=8)
		ax17.set_xlabel('Position',size=6)
		ax17.set_yticks(ax17.get_yticks()[::2])
		ax17.legend(loc=0,fontsize=5,labelspacing=0.05)

# 		Plot SD for MultiNucleotide Sequences
		ax18 = plt.subplot(gs[1])
		for dfNuc,lNuc in zip(ranMultiDataFrames,MultiNamesVal):
			elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()
			ax18.hist(elRegion,20,linewidth=0.3,label=lNuc,alpha=0.5)
		ax18.set_yticks(ax18.get_yticks()[::2])
		ax18.set_ylabel('Frequency',size=8)
		ax18.legend(loc=0,fontsize=5,labelspacing=0.1)
		ax18.set_xlabel('Standard Deviation Value',size=8)
		ax18.set_title('Standard Deviation for Polynucleotides within Random Region',size=8)

		sns.despine()
		plt.savefig(pp, format='pdf')
		print 'Plotted multinucleotide means'

# 	if set(['C','G','CG']) <= set(names):#set(['CG','C','G']).issubset(names):
# 		CGgroup,CGmean,CGstd = collect_sum_two_nucleotides(dfWindow,names,'C','G')
# 		ranCGgroup,ranCGmean,ranCGstd = collect_sum_two_nucleotides(ranWindow,names,'C','G')
# 
# 		gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
# 		gs.update(hspace=.8) # setting the space between the graphs
# 
# 		# Separate out those with only a double nucleotide search
# 		CpGNames = [names.index(i) for i in names if i == 'CG']
# 		CpGNamesVal = [names[i] for i in CpGNames]
# 		CpGDataFrames = [dfWindow[i] for i in CpGNames]
# 		ranCpGDataFrames = [ranWindow[i] for i in CpGNames]
# 
# 		# Plot actual C Gs available, to CpG presence
# 		ax19 = plt.subplot(gs[0],sharex=ax0)
# 		for dfNuc,lNuc in zip(CpGDataFrames,CpGNamesVal):
# 			ax19.plot(GlobalVariables.fillX,100*dfNuc.mean()/CGmean,linewidth=1,label='CG used for CpG in Element')
# 		for dfNuc,lNuc in zip(ranCpGDataFrames,CpGNamesVal):
# 			ax19.plot(GlobalVariables.fillX,dfNuc.mean()/ranCGmean*100,linewidth=1,label='CG used for CpG in Random Region')
# 		ax19.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
# 		ax19.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
# 		ax19.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
# 		ax19.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
# 		ax19.set_ylabel('% CG used for CpG',size=8)
# 		ax19.set_xlabel('Position',size=6)
# 		ax19.legend(loc=0,fontsize=5,labelspacing=0.1)
# 		ax19.set_title('% Available CG used for CpG',size=8)
# 		ax19.set_yticks(ax19.get_yticks()[::2])
# 		plt.xlim(0,GlobalVariables.num)
# 
# 		ax20 = plt.subplot(gs[1])
# 		for dfNuc,lNuc in zip(CpGDataFrames,CpGNamesVal):
# 			elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()#/CGstd*100
# 			ax20.hist(elRegion,20,linewidth=0.3,label='CG used for CpG in Element',alpha=0.5)
# 		for dfNuc,lNuc in zip(ranCpGDataFrames,CpGNamesVal):
# 			elRegion = dfNuc.loc[:,GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour].std()#/ranCGstd*100
# 			ax20.hist(elRegion,20,linewidth=0.3,label='CG used for CpG in Random Region',alpha=0.5)
# 		ax20.set_yticks(ax20.get_yticks()[::2])
# 		ax20.set_ylabel('Frequency',size=8)
# 		ax20.legend(loc=0,fontsize=5,labelspacing=0.1)
# 		ax20.set_xlabel('Standard Deviation Value',size=8)
# 		ax20.set_title('Standard Deviation within the Element',size=8)
# 		
# 		print 'Plotted the mean CG content and standard deviation'
	
	sns.despine()
	pp.savefig()
	pp.close()

def main(dfWindow,names,ranWindow,fileName):
	print 'Running GraphFangLibrary'
	graph_element_line_means(dfWindow,names,ranWindow,fileName)

if __name__ == "__main__":
	main()
