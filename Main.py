"""
Script to call and execute different sets of analyses for Element Structure

Wren Saylor
July 5 2017

To Do:
k means / scikit learn
dendrogram color clusters (what clustering parameters are used in literature, how do uces compare to other clustered groups)
cluster by difference between boundaries
snps, nucleosomes, replication origins/forks, tss
which are consitently = 
unittest

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
import ElementLibrary
import DirectionLibrary
import FangsLibrary
import MethylationLibrary
import RevCompLibrary
import TypeLibrary
import OverlapLibrary
import BinLibrary
import GraphFangLibrary
import GraphMethLibrary
import GraphSignalLibrary
import BokehLibrary
import GraphTableLibrary
import GraphClusterLibrary
import GraphDendrogramLibrary
import GraphKMeansLibrary
import pandas as pd
import RevCompOverlapLibrary
import BinLibrary
import GlobalVariables
import GraphMethExtendedLibrary

# set command line arguments
def get_args():
	# File lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=argparse.FileType('rU'),help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("mfile", type=argparse.FileType('rU'),help="A file containing a list of paths to the methylation files with unique names separated by newlines, with data for methylation position (chr, start,stop) and methylation % as fourth column'")
	parser.add_argument("rfile",type=str,help="A file containing the random regions equable with your elements to plot in contrast")

	# Genome Files
	parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
	parser.add_argument("-fa", "--fasta", type=str, default="hg19.fa")
	parser.add_argument("-o", "--overlapingelements", type=str, default="hg19_0based_exons.bed",help='A file containing all elements to check for overlaps in format "chr start stop", using for exons, but can be anything')#direction

	# Integer Parameters
	parser.add_argument("-t", "--total", type=int, default="600", help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e", "--element", type=int, default="200", help='size of your element (region without flanks), should be an even number')
	parser.add_argument("-i", "--inset", type=int, default="50", help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-w", "--window", type=int, default="11", help='size of sliding window, should be an odd number, previous studies have used 11')
	parser.add_argument("-b", "--bin", type=int, default="30", help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-mc", "--thresholdcoverage", type=int, default="10", help='size to threshold uncapped coverage of methylation data to send to % methylation, field often uses 10')
	parser.add_argument("-mp", "--thresholdpercentage", type=int, default="10", help='size to threshold % methylation data')
	parser.add_argument("-mf", "--methylationflank", type=int, default="10", help='The number of base pairs to look at outside of the element for the methylation clusterplots')

	# Specify which groups and graphs to run
	parser.add_argument('-type', "--elementype", default=['all'], nargs='*', choices=['all','intronic','exonic','intergenic'],help='which group types of element to run')
	parser.add_argument('-dir', "--elementdirection", default=[], nargs='*', choices=['+','-','='], help='which group direction of element to run')
	parser.add_argument('-rc', "--reversecomplement",action='store_true', help='if reverse complement sorting required')
	parser.add_argument('-p',"--plots",default=[],nargs='*',choices=['fang','methylation','signal','interactive','cluster','dendrogram','kmean','methextend'],help='the available graphs to plot')
	parser.add_argument('-nuc',"--nucleotideline",default=['A','T','C','G'],nargs='+',help='type the nucleotide string combinations to search for in the element')
	parser.add_argument('-str',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-align', "--elementalign",action='store_true', help='if want to align by exonic/intronic crossover')

	# Add additional descriptive file name information
	return parser.parse_args()

# for type and direction, separate the groups and run the analyses
def separate_dataframe_by_group(List,directionFeatures,typecolumn,fileName):
	# subset by bool presence
	bool = (directionFeatures[directionFeatures[typecolumn] == List])
	# if there is nothing in that set, skip
	if len(bool.index) != 0:
		Meth,dfWindow,Names = TypeLibrary.main(bool,fileName)
	return bool,Meth,dfWindow,Names

# the plotting options, if in the list of plot flags, run graph
def plot_graph_by_arg(pdMeth,rnMeth,dfWindow,names,ranWindow,fileName):
	if 'fang' in GlobalVariables.graphs:
		GraphFangLibrary.main(dfWindow,names,ranWindow,fileName)
	if 'signal' in GlobalVariables.graphs:
		GraphSignalLibrary.main(dfWindow,names,ranWindow,fileName)
	if 'methylation' in GlobalVariables.graphs:
		GraphMethLibrary.main(pdMeth,rnMeth,fileName)
	if 'methextend' in GlobalVariables.graphs:
		GraphMethExtendedLibrary.main(pdMeth,rnMeth,fileName)
	if 'interactive' in GlobalVariables.graphs:
		BokehLibrary.main(dfWindow,ranWindow,fileName)
	if 'cluster' in GlobalVariables.graphs:
		GraphClusterLibrary.main(dfWindow,ranWindow,pdMeth,rnMeth,names,fileName)
	if 'dendrogram' in GlobalVariables.graphs:
		GraphDendrogramLibrary.main(dfWindow,ranWindow,names,fileName)
	if 'kmean' in GlobalVariables.graphs:
		GraphKMeansLibrary.main(dfWindow,ranWindow,names,fileName)

# make a table for element and random region
def plot_chi_square_table(TableData,Title,ranTableData,ranTitle,fileName):
	GraphTableLibrary.main(TableData,Title,ranTableData,ranTitle,fileName)

# collect the counts for how many of each boundary direction for each type
def collect_counts_for_element_type(rangeFeatures):
	grouptype = rangeFeatures.groupby('type')['compareBoundaries'].value_counts()
	pdgroup = pd.DataFrame(grouptype)
	pdgroup.columns = ['counts']
	pvtable = pd.pivot_table(pdgroup,index=['type'],columns=['compareBoundaries'],values=['counts'])
	pvtable.columns.name = None
	pvtable.index.name = None
	pvtable.loc['all']= pvtable.sum()
	pvout = pvtable.xs('counts', axis=1, drop_level=True)
	return pvout

def main():
	# Collect arguments
	args = get_args()
	
	# Set global variables
	GlobalVariables.main(args)

	# for each element file provided
	for fileName in GlobalVariables.eFiles:
		# get the region info to work with
		rangeFeatures = ElementLibrary.main(fileName)
		# separate by direction
		directionFeatures,directionBins = DirectionLibrary.main(rangeFeatures,fileName)
		
		# for each random file provided
		randomFile = GlobalVariables.rFiles
		# Collect out filename labels 
		paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(GlobalVariables.uce,GlobalVariables.inuce,GlobalVariables.num,GlobalVariables.binDir,GlobalVariables.window,fileName,randomFile,GlobalVariables.stringName)
		
		# get the region info to work with
		randomFeatures = ElementLibrary.main(randomFile)
		# separate by direction
		randirFeatures,randirBins = DirectionLibrary.main(randomFeatures,randomFile)
		
		# Plot boundary probabilities
		BinLibrary.main(directionBins,randirBins,paramlabels)

		# Make tazble for the count of each direction for each type element
		elementGroup = collect_counts_for_element_type(rangeFeatures)
		randomGroup = collect_counts_for_element_type(randirFeatures)
		plot_chi_square_table(elementGroup,'Elements',randomGroup,'Random Regions','Direction_Count_{0}'.format(paramlabels))
		
		# All elements
		if 'all' in GlobalVariables.typeList:
			GlobalVariables.typeList.remove('all')
			pdMeth,allWindow,allNames = TypeLibrary.main(rangeFeatures,fileName)
			rnMeth,ranWindow,ranNames = TypeLibrary.main(randomFeatures,fileName)
			plot_graph_by_arg(pdMeth,rnMeth,allWindow,allNames,ranWindow,'all_{0}'.format(paramlabels))
			if GlobalVariables.revCom:
				revMeth,revWindow,revNames = RevCompLibrary.main(directionFeatures)
				ranrevMeth,ranrevWindow,ranrevNames = RevCompLibrary.main(randirFeatures)
				plot_graph_by_arg(revMeth,ranrevMeth,revWindow,revNames,ranrevWindow,'revComp_all_{0}'.format(paramlabels))

# 		# By Type
		for type in GlobalVariables.typeList:
			typeBool,typeMeth,typeWindow,typeNames = separate_dataframe_by_group(type,directionFeatures,'type',fileName)
			rantypeBool,rantypeMeth,rantypeWindow,rantypeNames = separate_dataframe_by_group(type,randirFeatures,'type',randomFile)
			plot_graph_by_arg(typeMeth,rantypeMeth,typeWindow,typeNames,rantypeWindow,'{0}_{1}'.format(type,paramlabels))
			if GlobalVariables.revCom:
				typercMeth,typercWindow,typercNames = RevCompLibrary.main(typeBool)
				rantypercMeth,rantypercWindow,rantypercNames = RevCompLibrary.main(rantypeBool)
				plot_graph_by_arg(typercMeth,rantypercMeth,typercWindow,typercNames,rantypercWindow,'revComp_{0}_{1}'.format(type,paramlabels))

		# By Direction
		for dir in GlobalVariables.dirList:
			dirBool,dirMeth,dirWindow,dirNames = separate_dataframe_by_group(dir,directionFeatures,'compareBoundaries',fileName)
			randirBool,randirMeth,randirWindow,randirNames = separate_dataframe_by_group(dir,randirFeatures,'compareBoundaries',randomFile)
			plot_graph_by_arg(dirMeth,randirMeth,dirWindow,dirNames,randirWindow,'all_{0}_{1}'.format(dir,paramlabels))
			for type in GlobalVariables.typeList:
				typeBool,typeMeth,typeWindow,typeNames = separate_dataframe_by_group(type,directionFeatures,'type',fileName)
				rantypeBool,rantypeMeth,rantypeWindow,rantypeNames = separate_dataframe_by_group(type,randirFeatures,'type',randomFile)
				plot_graph_by_arg(typeMeth,rantypeMeth,typeWindow,typeNames,rantypeWindow,'{0}_{1}_{2}'.format(type,dir,paramlabels))

		# Re-align by exon/intron crossover
		if GlobalVariables.overlapInset:
			ranstartcrossboundary,ranendcrossboundary,rancompleteelement,raninteriorelement,ranoverlapTable = OverlapLibrary.main(randomFeatures,fileName)
			rancrossMeth,rancrossWindow,rancrossNames = RevCompOverlapLibrary.main(ranstartcrossboundary,ranendcrossboundary)
			startcrossboundary,endcrossboundary,completeelement,interiorelement,overlapTable = OverlapLibrary.main(rangeFeatures,fileName)
			crossMeth,crossWindow,crossNames  = RevCompOverlapLibrary.main(startcrossboundary,endcrossboundary)
			completeelement.name = 'WithinExon'
			interiorelement.name = 'ContainsExon'
			plot_chi_square_table(overlapTable.T,'Exonic Elements',ranoverlapTable.T,'Exonic Random Regions','Exonic_overlap_Count_{0}'.format(paramlabels))
			plot_graph_by_arg(crossMeth,rancrossMeth,crossWindow,crossNames,rancrossWindow,'align_Crossboundary_{0}'.format(paramlabels))
			for element,random in zip([completeelement,interiorelement],[rancompleteelement,raninteriorelement]):
				alignMeth,alignWindow,alignNames = TypeLibrary.main(element,fileName)
				ranalignMeth,ranalignWindow,ranalignNames = TypeLibrary.main(random,fileName)
				plot_graph_by_arg(alignMeth,ranalignMeth,alignWindow,alignNames,ranalignWindow,'align_{0}_{1}'.format(element.name,paramlabels))

if __name__ == "__main__":
	main()
