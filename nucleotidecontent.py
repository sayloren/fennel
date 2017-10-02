"""
Script to run mean AT conent analysis

Wren Saylor
September 1 2017

To Do:
unittest
change directionality arg column
streamline and remove excess columns
.index len to break up large datasets

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
import pandas as pd
from collections import defaultdict
import pybedtools as pbt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import seaborn as sns
import math
import matplotlib as mpl
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
from scipy.stats import chisquare
from scipy.interpolate import splrep, splev
import scipy.stats as ss
from scipy.stats import mstats

# set command line arguments
def get_args():
	# File lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=str,help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("-r","--randomfile",type=argparse.FileType('rU'),help="A file containing a list of paths to the random regions equable with your elements to plot in contrast") # option to not use random at all, or generate randoms?

	# Columns in element file - all 0 based
	parser.add_argument("-lc", "--labelcolumn",type=int,help='column in the element file where the label (exonic, intergenic, intronic) is')
	parser.add_argument("-dc", "--directionalitycolumn",type=int,help='column in the element file where directionality is, if not supplied will infer by AT content')
	parser.add_argument("-ic", "--idcolumn",type=int,help='column in the element file where the id is, if not provided will be generated for the sliding window collection')

	# Genome Files
	parser.add_argument("-g","--genome",type=str, default="hg19.genome")
	parser.add_argument("-f","--fasta",type=str,default="hg19.fa")

	# Integer Parameters
	parser.add_argument("-t","--total",type=int,default="600",help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e","--element",type=int,help='size of your element (region without flanks), should be an even number, if not provided will use the smallest size of your input data')# ,default="200"
	parser.add_argument("-i","--inset",type=int,default="50",help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-w","--window",type=int,default="11",help='size of sliding window, should be an odd number, previous studies have used 11')
	parser.add_argument("-b","--bin",type=int,default="100",help='size of bins used to compare element ends and determine directionality')

	# Plot parameters
	parser.add_argument('-s',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-c',"--reversecomplement",action='store_true',help='if you want the reverse complement to be plotted')
	parser.add_argument("-n", "--numberrandomassignments",type=int,default="1",help='the number of times to to randomly assign direction, will only be used when "--reversecomplement" is "random"')

	# Add additional descriptive file name information
	return parser.parse_args()

# set all args that will be used throughout the script
def set_global_variables(args):
	# Integer parameters
	global num
	global elementsize
	global inuce
	global window
	global binDir
	global halfwindow
	global fillX
	num = args.total
	elementsize = args.element
	inuce = args.inset
	window = args.window
	binDir = args.bin
	halfwindow = ((window/2)+1)
	fillX = range(0,(num-window))

	# Element, random regions and methylation files
	global eFiles
	global rFiles
	eFiles = args.efile
	if args.randomfile:
		rFiles = [line.strip() for line in args.randomfile]
	else:
		rFiles = None
	
	# Column labels in element file
	global labelcolumn
	global directionalitycolumn
	global idcolumn
	labelcolumn = args.labelcolumn
	directionalitycolumn = args.directionalitycolumn
	idcolumn = args.idcolumn

	# Genome files from UCSC
	global sizeGenome
	global faGenome
	sizeGenome = args.genome
	faGenome = args.fasta

	# Lists with the types and directions to use
	global nucList
	nucList = ['A','T']
	
	# A string to add to the out file name in case you want to set up runs and let be
	global stringName
	stringName = args.stringname
	
	global reverseComplement
	global randomassignments
	reverseComplement = args.reversecomplement
	randomassignments = args.numberrandomassignments

def set_ploting_parameters():
	# Locations for plotting with sliding window
	global plotLineLocationOne # upstream element boundary
	global plotLineLocationTwo # downstream element boundary
	global plotLineLocationThree # upstream element inset
	global plotLineLocationFour # downstream element inset
	plotLineLocationOne = (((num-uce)/2)+(inuce-halfwindow))
	plotLineLocationTwo = (((num-uce)/2)+(uce-inuce-halfwindow))
	plotLineLocationThree = (((num-uce)/2)-halfwindow)
	plotLineLocationFour = (((num-uce)/2)+uce-halfwindow)
	
	# Locations for plotting without the sliding window
	global plotLineLocationOneFull # upstream element boundary
	global plotLineLocationTwoFull # downstream element boundary
	global plotLineLocationThreeFull # upstream element inset
	global plotLineLocationFourFull # downstream element inset
	plotLineLocationOneFull = (((num-uce)/2)+inuce)
	plotLineLocationTwoFull = (((num-uce)/2)+(uce-inuce))
	plotLineLocationThreeFull = ((num-uce)/2)
	plotLineLocationFourFull = (((num-uce)/2)+uce)

# get bt features
def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# get the correct range for fang evaluation
def collect_coordinates_for_element_positions(btFeatures):
	midFeatures = pd.read_table(btFeatures.fn, header=None)
	midFeatures['middle'] = midFeatures.loc[:,1:2].mean(axis=1).astype(int)
	midFeatures['chr'] = midFeatures.loc[:,0].astype(str)
	midFeatures['start'] = midFeatures.loc[:,1].astype(int)
	midFeatures['end'] = midFeatures.loc[:,2] .astype(int)
	midFeatures['size'] = midFeatures.loc[:,2].astype(int)-midFeatures.loc[:,1].astype(int)
	global uce
	if elementsize:
		uce = elementsize
	else:
		getmin = midFeatures['size'].min()
		if (getmin % 2) == 0:
			uce = getmin
		else:
			uce = getmin + 1
	flankSize = (num - uce)/2
	inregion = uce-(inuce*2)
	midFeatures['sCenter'] = midFeatures['middle'].astype(int) - (inregion/2)
	midFeatures['eCenter'] = midFeatures['middle'].astype(int) + (inregion/2)
	midFeatures['sEdge'] = midFeatures.loc[:,1].astype(int) + inuce
	midFeatures['eEdge'] = midFeatures.loc[:,2].astype(int) - inuce
	midFeatures['sBoundary'] = midFeatures.loc[:,1].astype(int) - flankSize
	midFeatures['eBoundary'] = midFeatures.loc[:,2].astype(int) + flankSize
	if idcolumn:
		midFeatures['id'] = midFeatures.loc[:,idcolumn]
	else:
		midFeatures.insert(0,'id',range(0,0+len(midFeatures)))
	if labelcolumn:
		midFeatures['type'] = midFeatures.loc[:,labelcolumn]
	if directionalitycolumn:
		midFeatures['directionality'] = midFeatures.loc[:,directionalitycolumn]
	return midFeatures

# check that the coordinates with the surrounding regions don't fall off the end of genome, if they do, print and skip
def check_coords_beyond_genome(rangeFeatures):
	genomeBedtool = get_bedtools_features(sizeGenome)
	genomeFeatures = pd.read_table(genomeBedtool.fn, header=None)
	genomeFeatures['chr'] = genomeFeatures.loc[:,0]
	genomeFeatures['end'] = genomeFeatures.loc[:,1]
	initiallength=len(rangeFeatures.index)
	chrList = rangeFeatures['chr'].unique()
	for chr in chrList:
		end = genomeFeatures.loc[(genomeFeatures['chr'] == chr),'end'].values[0]
		rangeFeatures = rangeFeatures[(rangeFeatures['chr'] == chr) & (rangeFeatures['eBoundary'] < end)]
		# can check if any starts are negative
	checklength = initiallength - len(rangeFeatures.index)
	print "There were {0} features where the surrounding features streached beyond the end genome boundary".format(checklength)
	return rangeFeatures

# get the strings for sliding window regions
def get_fasta_for_element_coordinates(rangeFeatures):
# 	rangeFeatures = check_coords_beyond_genome(rangeFeatures)
	rangeFeatures['sBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sBoundary','start']].values.astype(str).tolist()))
	rangeFeatures['sEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','sEdge']].values.astype(str).tolist()))
	rangeFeatures['MiddleSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sCenter','eCenter']].values.astype(str).tolist()))
	rangeFeatures['eEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','eEdge','end']].values.astype(str).tolist()))
	rangeFeatures['eBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','end','eBoundary']].values.astype(str).tolist()))
	rangeFeatures['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','end']].values.astype(str).tolist()))
	rangeFeatures['combineString'] = rangeFeatures['sBoundarySeq'].astype(str) + rangeFeatures['sEdgeSeq'].astype(str) + rangeFeatures['MiddleSeq'].astype(str) + rangeFeatures['eEdgeSeq'].astype(str) + rangeFeatures['eBoundarySeq'].astype(str)
	rangeFeatures['combineString'] = rangeFeatures['combineString'].str.upper()
	rangeFeatures['feature'] = rangeFeatures['feature'].str.upper()
	rangeFeatures['reverseComplement'] = rangeFeatures.apply(lambda row: reverse_complement_dictionary(row['combineString']),axis=1)
	rangeFeatures.drop(['sBoundarySeq','sEdgeSeq','MiddleSeq','eEdgeSeq','eBoundarySeq','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary'],axis=1)
	return rangeFeatures

# get the reverse complement
def reverse_complement_dictionary(sequence):
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	return "".join([seqDict[base] for base in reversed(sequence)])

# used in get_fasta_for_element_coordinates to extract just the fasta strings
def get_just_fasta_sequence_for_feature(inFeature):
	seqFeature = inFeature.sequence(fi=faGenome)
	outFeature = pd.read_table(seqFeature.seqfn)
	outSequence = outFeature[::2]
	outSequence = outSequence.reset_index(drop=True)
	return outSequence

# Get the percentage AT in the element
def percentage_at_for_element(region,group):
	collectAT = []
	for r in region:
		collectAT.append(eval('100*float(r.count("A") + r.count("a") + r.count("T") + r.count("t"))/len(r)'))
	pdAT = pd.DataFrame(collectAT)
	print 'Mean AT content for {0} elements in {1} is {2} Percent'.format(len(region.index),group,pdAT.mean())

# get coordinates with flanking regions
def collect_element_coordinates(fileName):
	btFeatures = get_bedtools_features(fileName)
	subsetFeatures = collect_coordinates_for_element_positions(btFeatures)
	rangeFeatures = get_fasta_for_element_coordinates(subsetFeatures)
	return rangeFeatures

# Do the comparison between boundaries to get + - or =
def calculate_nucleotides_at(element,size):
	start = element[:size]
	end = element[-size:]
	perSize = []
	perSize.append(eval('100*float(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
	perSize.append(eval('100*float(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
	return perSize

# Directionality, as inferred by comparing first and last n base pairs from input parameters
def compare_boundaries_size_n(element,size):
	perSize = calculate_nucleotides_at(element,size)
	# give + - = depending on which side has larger AT content
	if perSize[0] > perSize[1]: outList = '+'
	if perSize[1] > perSize[0]: outList = '-'
	if perSize[1] == perSize[0]: outList = '='
	return outList

# With the results from compare_boundaries_size_n per each element, evaluate directionality into new column
def evaluate_boundaries_size_n(rangeFeatures,fileName):
	rangeFeatures['compareBoundaries'] = rangeFeatures.apply(lambda row: (compare_boundaries_size_n(row['feature'],binDir)),axis=1)
	compareEnds = pd.DataFrame(rangeFeatures[['chr','start','end','compareBoundaries']])
	print 'Sorting the element boundaries by bin size {0}'.format(binDir)
	return rangeFeatures

# run the sliding window for each nucleotide string
def run_sliding_window_for_each_nucleotide_string(features,label):
	outCollect = []
	for element,id in zip(features,label):
		outElement = {id: []}
		outList = {key:[] for key in nucList}
		n = num
		s = 1 # size to jump for sliding window
		start, end = 0, window
		while end < n:
			current = element[start:end]
			for key in nucList:
				percentage = float(100*current.count(key)/len(current))
				outList[key].append(percentage)
			start, end = start + s, end + s
		outElement[id].append(outList)
		outCollect.append(outElement)
	return outCollect

# convert to a data frame with a list for each element x nucleotide string
def flatten_data_from_sliding_window(outCollect):
	outFlatten = pd.DataFrame()
	outIndex = []
	for d in outCollect:
		for key,values in d.items():
			outFlatten = outFlatten.append(values,ignore_index=True)
			outIndex.append(key)
	outFlatten.index = outIndex
	return outFlatten

# turn each list of element x nucleotide string into a separate df, within a larger df
def convert_sliding_window_to_dataframe(outFlatten):
	names = outFlatten.columns.tolist()
	collectNucDF = []
	for nuc in names:
		nuc = outFlatten[[nuc]]
		nuc.columns = ['temp']
		split = pd.DataFrame(nuc['temp'].values.tolist(),index=nuc.index)
		collectNucDF.append(split)
	return collectNucDF,names

# Put all the random regions for each string search into one df, for stat collection
def sliding_window_df_to_collect_all_random(collectRandom,allNames):
	dictNames = {key:[] for key in allNames}
	for random in collectRandom:
		for df,key in zip(random,dictNames):
			dictNames[key].append(df)
	catCollect = []
# 	dictCollect = {key:[] for key in allNames}
	for key,values in dictNames.items():
		cat = pd.concat(values)
		catCollect.append(cat)
# 		dictCollect[key].append(cat)
	return catCollect

# Wrapper for running sliding window and converting it into easy to use format
def sliding_window_wrapper(features,label):
	outCollect = run_sliding_window_for_each_nucleotide_string(features,label)
	# outCollect is a list containing a dictionary for each nucleotide string search, with a list of the % at each search position
	outFlatten=flatten_data_from_sliding_window(outCollect)
	# outFlatten is a data frame with a column for each search string, containing a list of % at each search position
	outDataFrame,names = convert_sliding_window_to_dataframe(outFlatten)
	# outDataFrame is a list of dataframes for each search string, with the each % at each search position in a separate column
	print 'Retrieved sliding window data for nucleotides strings {0}'.format(names)
	return outDataFrame,names

# Get group, mean and standard deviation for AT
def collect_sum_two_nucleotides(dfWindow,names,nuc1,nuc2):
	ATNames = [names.index(i) for i in names if nuc1 in i or nuc2 in i]
	ATDataFrames = [dfWindow[i] for i in ATNames]
	ATconcat = pd.concat(ATDataFrames,axis=1)
	ATgroup = ATconcat.groupby(ATconcat.columns,axis=1).sum()
	ATmean = ATgroup.mean()
	ATstd = ATgroup.std()
	return ATgroup, ATmean, ATstd

# Make graphs for fangs, with rc sorting
def graph_element_line_means_with_rc_sorted(dfWindow,names,revWindow,fileName,collectRandom,collectRandomRC,denseRandom,denseRandomRC):
	set_ploting_parameters()
	ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(denseRandom,names,'A','T')
	revATgroup,revATmean,revATstd = collect_sum_two_nucleotides(revWindow,names,'A','T')
	revranATgroup,revranATmean,revranATstd = collect_sum_two_nucleotides(denseRandomRC,names,'A','T')
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"
	sns.set_style('ticks')
	gs = gridspec.GridSpec(2,2,height_ratios=[3,1])
	gs.update(hspace=.8)
	pp = PdfPages('Fangs_{0}.pdf'.format(fileName))
	plt.figure(figsize=(14,7))
	plt.suptitle(info,fontsize=10)
	sns.set_palette("husl",n_colors=8)

	# Stats
# 	wilcoxonsignedrank = ss.wilcoxon(ATmean,ranATmean)
# 	wilcoxonsignedrankreverse = ss.wilcoxon(revATmean,revranATmean)
# 	sdelement = ATgroup.loc[:,plotLineLocationThree:plotLineLocationFour].std()
# 	sdrandom = ranATgroup.loc[:,plotLineLocationThree:plotLineLocationFour].std()
# 	sdelementreverse = revATgroup.loc[:,plotLineLocationThree:plotLineLocationFour].std()
# 	sdrandomreverse = revranATgroup.loc[:,plotLineLocationThree:plotLineLocationFour].std()
# 	sdkruskal = mstats.kruskalwallis(sdelement,sdrandom)
# 	sdkruskalreverse = mstats.kruskalwallis(sdelementreverse,sdrandomreverse)
# 	ax0.text(20,90,'Wilcox Signed Rank P-value {:0.1e}'.format(wilcoxPSRMean[1]),size=6,clip_on=False)
# 	ax2.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)

	ax0 = plt.subplot(gs[0,0])
	for dfNuc in collectRandom:
		ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(dfNuc,names,'A','T')
		ax0.plot(fillX,ranATmean,linewidth=1,alpha=0.3)
	ax0.plot(fillX,ATmean,linewidth=1,label='AT element')
# 	ax0.fill_between(fillX,ATmean+ATstd,ATmean-ATstd,label='',alpha=0.2)
	ax0.axvline(x=plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.axvline(x=plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
# 	ax0.hlines(y=66,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
# 	ax0.text(32,65,'11bp sliding window',size=6)
	ax0.set_ylabel('% AT Content',size=8)
	ax0.set_xlabel('Position',size=6)
	ax0.set_title('Mean AT Content With Standard Deviation',size=8)
	ax0.set_yticks(ax0.get_yticks()[::2])
	plt.xlim(0,num)

	ax1 = plt.subplot(gs[1,0],sharex=ax0)
	for dfNuc in collectRandom:
		ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(dfNuc,names,'A','T')
		ax1.plot(fillX,ranATstd,linewidth=1,alpha=0.3)
	ax1.plot(fillX,ATstd,linewidth=1,label='AT element')
	ax1.axvline(x=plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.axvline(x=plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.set_yticks(ax1.get_yticks()[::2])
	ax1.set_xlabel('Position',size=6)
	ax1.set_ylabel('SD',size=8)
	ax1.set_title('Standard Deviation',size=8)
	plt.setp(ax1.get_xticklabels(), visible=True)
	
	sns.set_palette("husl",n_colors=8)
	ax2 = plt.subplot(gs[0,1])
	for rcNuc in collectRandomRC:
		ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(rcNuc,names,'A','T')
		ax2.plot(fillX,ranATmean,linewidth=1,alpha=0.3)
	ax2.plot(fillX,revATmean,linewidth=1,label='AT element')
# 	ax0.fill_between(fillX,revATmean+revATstd,revATmean-revATstd,label='',alpha=0.2)
	ax2.axvline(x=plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax2.axvline(x=plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax2.axvline(x=plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax2.axvline(x=plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
# 	ax2.hlines(y=66,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
# 	ax2.text(32,65,'11bp sliding window',size=6)
	ax2.set_ylabel('% AT Content',size=8)
	ax2.set_xlabel('Position',size=6)
	ax2.set_title('Mean AT Content With Standard Deviation',size=8)
	ax2.set_yticks(ax2.get_yticks()[::2])
	plt.xlim(0,num)

	ax3 = plt.subplot(gs[1,1],sharex=ax0)
	for rcNuc in collectRandomRC:
		ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(rcNuc,names,'A','T')
		ax3.plot(fillX,ranATstd,linewidth=1,alpha=0.3)
	ax3.plot(fillX,revATstd,linewidth=1,label='AT element')
	ax3.axvline(x=plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax3.axvline(x=plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax3.axvline(x=plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax3.axvline(x=plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax3.set_yticks(ax1.get_yticks()[::2])
	ax3.set_xlabel('Position',size=6)
	ax3.set_ylabel('SD',size=8)
	ax3.set_title('Standard Deviation',size=8)
	plt.setp(ax1.get_xticklabels(), visible=True)

	sns.despine()
	pp.savefig()
	pp.close()

# Make graphs for fangs
def graph_element_line_means(dfWindow,names,fileName,Random,denseRandom):
	set_ploting_parameters()
	ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(denseRandom,names,'A','T')
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"
	sns.set_style('ticks')
	gs = gridspec.GridSpec(2,1,height_ratios=[3,1])
	gs.update(hspace=.8)
	pp = PdfPages('Fangs_{0}.pdf'.format(fileName))
	plt.figure(figsize=(14,7))
	plt.suptitle(info,fontsize=10)
	sns.set_palette("husl",n_colors=8)
	
	# Stats
# 	wilcoxonsignedrank = ss.wilcoxon(ATmean,ranATmean)
# 	sdelement = ATgroup.loc[:,plotLineLocationThree:plotLineLocationFour].std()
# 	sdrandom = ranATgroup.loc[:,plotLineLocationThree:plotLineLocationFour].std()
# 	sdkruskal = mstats.kruskalwallis(sdelement,sdrandom)
# 	ax0.text(20,90,'Wilcox Signed Rank P-value {:0.1e}'.format(wilcoxPSRMean[1]),size=6,clip_on=False)
# 	ax2.text(16.25,14.5,'KW P-value {:0.1e}'.format(kruskalSD[1]),size=6,clip_on=False)

	ax0 = plt.subplot(gs[0,:])
	ax1 = plt.subplot(gs[1,:],sharex=ax0)
	for dfNuc in Random:
		ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(dfNuc,names,'A','T')
		ax0.plot(fillX,ranATmean,linewidth=1,alpha=0.3)
		ax1.plot(fillX,ranATstd,linewidth=1,alpha=0.3)
	ax0.plot(fillX,ATmean,linewidth=1,label='AT element')
	ax0.axvline(x=plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.axvline(x=plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
# 	ax0.hlines(y=66,xmin=20,xmax=31,linewidth=.5,color='#081d58',zorder=0)
# 	ax0.text(32,65,'11bp sliding window',size=6)
	ax0.set_ylabel('% AT Content',size=8)
	ax0.set_xlabel('Position',size=6)
	ax0.set_title('Mean AT Content With Standard Deviation',size=8)
	ax0.set_yticks(ax0.get_yticks()[::2])
	plt.xlim(0,num)
	ax1.plot(fillX,ATstd,linewidth=1,label='AT element')
	ax1.axvline(x=plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.axvline(x=plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.set_yticks(ax1.get_yticks()[::2])
	ax1.set_xlabel('Position',size=6)
	ax1.set_ylabel('SD',size=8)
	ax1.set_title('Standard Deviation',size=8)
	plt.setp(ax1.get_xticklabels(),visible=True)

	sns.despine()
	pp.savefig()
	pp.close()

# For type groups, separate the groups and run the analyses
def separate_dataframe_by_group(List,directionFeatures,typecolumn,fileName):
	# subset by bool presence
	bool = (directionFeatures[directionFeatures[typecolumn] == List])
	dfWindow,Names = sliding_window_wrapper(bool['combineString'],bool['id'])
	return bool,dfWindow,Names

# Sliding window RCsorting
def sort_sliding_window_by_directionality(negStr,posStr):
	negDF, negNames = sliding_window_wrapper(negStr['reverseComplement'],negStr['id'])
	posDF, posNames = sliding_window_wrapper(posStr['combineString'],posStr['id'])
	compWindow = []
	for x, y in zip(negDF,posDF):
		tempCat = pd.concat([x,y],axis=1)
		tempGroup = tempCat.groupby(tempCat.columns,axis=1).sum()
		compWindow.append(tempGroup)
	return compWindow, negNames

# Separate and sort by plus and minus orientation
def sort_elements_by_directionality(directionFeatures,columnCompare):
	negStr = (directionFeatures[(directionFeatures[columnCompare] == '-')])
	posStr = (directionFeatures[(directionFeatures[columnCompare] == '+')])
	compWindow, compNames = sort_sliding_window_by_directionality(negStr,posStr)
	return compWindow,compNames

# Get the empirical probability for each direction classification
def make_probabilites_for_direction(directionFeatures,probabilitycolumn):
	lenAll = float(len(directionFeatures.index))
	numPlus = (directionFeatures[probabilitycolumn] == '+').sum()/lenAll
	numMinus = (directionFeatures[probabilitycolumn] == '-').sum()/lenAll
	numEqual = (directionFeatures[probabilitycolumn] == '=').sum()/lenAll
	probOptions = [numMinus,numPlus,numEqual]
	return probOptions

def main():
	# Get and set parameters
	args = get_args()
	set_global_variables(args)
	paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(elementsize,inuce,num,binDir,window,eFiles,stringName)
	
	# Get coords and strings for elements
	rangeFeatures = collect_element_coordinates(eFiles)
	percentage_at_for_element(rangeFeatures['feature'],eFiles)
	directionFeatures = evaluate_boundaries_size_n(rangeFeatures,eFiles)
	
	# Get the probability for each directional assignment, and use to randomly assign the correct number of random directions
	dirOptions = ['-','+','=']
	if directionalitycolumn:
		probOptions = make_probabilites_for_direction(directionFeatures,'directionality')
	else:
		probOptions = make_probabilites_for_direction(directionFeatures,'compareBoundaries')
	
	if not labelcolumn:
		allWindow,allNames = sliding_window_wrapper(rangeFeatures['combineString'],rangeFeatures['id'])
		if reverseComplement:
			if directionalitycolumn:
				revWindow,revNames = sort_elements_by_directionality(directionFeatures,'directionality')
			else:
				revWindow,revNames = sort_elements_by_directionality(directionFeatures,'compareBoundaries')
		spreadRandom = []
		spreadRandomRC = []
		if rFiles:
			for randomFile in rFiles:
				randomFeatures = collect_element_coordinates(randomFile)
				randirFeatures= evaluate_boundaries_size_n(randomFeatures,randomFile)
				sWindow,sNames = sliding_window_wrapper(randirFeatures['combineString'],randirFeatures['id'])
				spreadRandom.append(sWindow)
				if reverseComplement:
					ranrevWindow, ranrevNames = sort_elements_by_directionality(randirFeatures,'compareBoundaries')
					spreadRandomRC.append(ranrevWindow)
		else:
			for i in range(randomassignments):
				directionFeatures['randomDirection'] = np.random.choice(dirOptions,len(directionFeatures.index),p=probOptions)
				randirWindow, randirNames = sort_elements_by_directionality(directionFeatures,'randomDirection')
				spreadRandom.append(randirWindow)
				spreadRandomRC.append(randirWindow)
		denseRandom = sliding_window_df_to_collect_all_random(spreadRandom,allNames)
		if reverseComplement:
			denseRandomRC = sliding_window_df_to_collect_all_random(spreadRandomRC,allNames)
			graph_element_line_means_with_rc_sorted(allWindow,allNames,revWindow,'all_rc_{0}'.format(paramlabels),spreadRandom,spreadRandomRC,denseRandom,denseRandomRC)
		else:
			graph_element_line_means(allWindow,allNames,'all_{0}'.format(paramlabels),spreadRandom,denseRandom)
	if labelcolumn:
		typeList = rangeFeatures[labelcolumn].unique()
		for type in typeList:
			typeBool,typeWindow,typeNames = separate_dataframe_by_group(type,rangeFeatures,'type',eFiles)
			percentage_at_for_element(typeBool['feature'],'{0}, {1}'.format(eFiles,type))
			if directionalitycolumn:
				probOptionstype = make_probabilites_for_direction(typeBool,'directionality')
			else:
				probOptionstype = make_probabilites_for_direction(typeBool,'compareBoundaries')
			if reverseComplement:
				if directionalitycolumn:
					typeWindowRC,typeNamesRC = sort_elements_by_directionality(typeBool,'directionality')
				else:
					typeWindowRC,typeNamesRC = sort_elements_by_directionality(typeBool,'compareBoundaries')
			spreadRandomtype = []
			spreadRandomtypeRC = []
			if rFiles:
				for randomFile in rFiles:
					randomFeatures = collect_element_coordinates(randomFile)
					randirFeatures= evaluate_boundaries_size_n(randomFeatures,randomFile)
					rantypeBool,rantypeWindow,rantypeNames = separate_dataframe_by_group(type,randirFeatures,'type',randomFile)
					spreadRandomtype.append(rantypeWindow)
					if reverseComplement:
						typeWindowRCrandom,typeNamesRCrandom = sort_elements_by_directionality(rantypeBool,'compareBoundaries')
						spreadRandomtypeRC.append(typeWindowRCrandom)
			else:
				for i in range(randomassignments):
					typeBool['randomDirectiontype'] = np.random.choice(dirOptions,len(typeBool.index),p=probOptionstype)
					typedirWindow, typedirNames = sort_elements_by_directionality(typeBool,'randomDirectiontype')
					spreadRandomtype.append(typedirWindow)
					spreadRandomtypeRC.append(typedirWindow)
			denseRandomtype = sliding_window_df_to_collect_all_random(spreadRandomtype,typeNames)
			if reverseComplement:
				denseRandomtypeRC = sliding_window_df_to_collect_all_random(spreadRandomtypeRC,allNames)
				graph_element_line_means_with_rc_sorted(typeWindow,typeNames,typeWindowRC,'{0}_rc_{1}'.format(type,paramlabels),spreadRandomtype,spreadRandomtypeRC,denseRandomtype,denseRandomtypeRC)
			else:
				graph_element_line_means(typeWindow,typeNames,'{0}_{1}'.format(type,paramlabels),spreadRandomtype,denseRandomtype)

if __name__ == "__main__":
	main()