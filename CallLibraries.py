"""
Script to call and execute different sets of analyses for Element Structure

Wren Saylor
July 5 2017

To Do:
k means / scikit learn
dendrogram color clusters (what clustering parameters are used in literature, how do uces compare to other clustered groups)
cluster by difference between boundaries
snps, nucleosomes, replication origins/forks, tss

more constrined = boundary graph from emperical spread
c/g available
AT balance random
exons - split intron/intergenic, include exon direcitonality, move cross boundary to tss
which are consitently = 
AT liklihood (emperical)

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

# set command line arguments
def get_args():
	# File lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=argparse.FileType('rU'), help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("mfile", type=argparse.FileType('rU'), help="A file containing a list of paths to the methylation files with unique names separated by newlines, with data for methylation position (chr, start,stop) and methylation % as fourth column'")
	parser.add_argument("rfile",type=argparse.FileType('rU'), help="A file containing a list of paths to the random regions equable with your elements to plot in contrast")

	# Genome Files
	parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
# 	parser.add_argument("-n", "--nucleosome", type=str, help="A bedgraph file with data for nucleosome positions, form 'chr, start, stop, occupancy'")
# 	parser.add_argument("-s", "--snp", type=str, help="A file with data for snps, form 'chr, start, stop(start+size alt-mutated), ref, ref_size, alt, alt_size, af_adj'")
	parser.add_argument("-fa", "--fasta", type=str, default="hg19.fa")
	parser.add_argument("-o", "--overlapingelements", type=str, default="hg19_0based_exons.bed",help='A file containing all elements to check for overlaps in format "chr start stop", using for exons, but can be anything')#direction

	# Integer Parameters
	parser.add_argument("-t", "--total", type=int, default="600", help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e", "--element", type=int, default="200", help='size of your element (region without flanks), should be an even number')
	parser.add_argument("-i", "--inset", type=int, default="50", help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-w", "--window", type=int, default="11", help='size of sliding window, should be an odd number, previous studies have used 11')
	parser.add_argument("-b", "--bin", type=int, default="30", help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-mc", "--thresholdcoverage", type=int, default="10", help='size to threshold uncapped coverage of methylation data to send to % methylation, field often uses 10')
	parser.add_argument("-mp", "--thresholdpercentage", type=int, default="0", help='size to threshold % methylation data')
	parser.add_argument("-mf", "--methylationflank", type=int, default="20", help='The number of base pairs to look at outside of the element for the methylation clusterplots')

	# Specify which groups and graphs to run
	parser.add_argument('-type', "--elementype", default=[], nargs='*', choices=['all','intronic','exonic','intergenic'],help='which group types of element to run')
	parser.add_argument('-dir', "--elementdirection", default=[], nargs='*', choices=['+','-','='], help='which group direction of element to run')
	parser.add_argument('-rc', "--reversecomplement",action='store_true', help='if reverse complement sorting required')
	parser.add_argument('-p',"--plots",default=[],nargs='*',choices=['fang','methylation','signal','interactive','cluster','dendrogram','kmean'],help='the available graphs to plot')
	parser.add_argument('-nuc',"--nucleotideline",default=['A','T'],nargs='+',help='type the nucleotide string combinations to search for in the element')
	parser.add_argument('-str',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-align', "--elementalign",action='store_true', help='if want to align by exonic/intronic crossover')

	# Add additional descriptive file name information
	return parser.parse_args()

# for type and direction, separate the groups and run the analyses
def groupSeparate(List,directionFeatures,typecolumn,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs,typecol):

	# subset by bool presence
	bool = (directionFeatures[directionFeatures[typecolumn] == List])

	# if there is nothing in that set, skip
	if len(bool.index) != 0:
		Meth,dfWindow,Names = TypeLibrary.main(bool,typecol,fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
	return bool,Meth,dfWindow,Names

# the plotting options, if in the list of plot flags, run graph
def plotGraphs(pdMeth,rnMeth,dfWindow,names,ranWindow,fileName,num,uce,inuce,window,graphs,nucLine,methFlank):
	if 'fang' in graphs:
		GraphFangLibrary.main(dfWindow,names,ranWindow,fileName,num,uce,inuce,window,nucLine)
	if 'signal' in graphs:
		GraphSignalLibrary.main(dfWindow,names,ranWindow,fileName,num,uce,inuce,window,nucLine)
	if 'methylation' in graphs:
		GraphMethLibrary.main(pdMeth,rnMeth,fileName,num,uce,inuce,window)
	if 'interactive' in graphs:
		BokehLibrary.main(dfWindow,ranWindow,fileName,num,uce,inuce,window,nucLine)
	if 'cluster' in graphs:
		GraphClusterLibrary.main(dfWindow,ranWindow,pdMeth,rnMeth,names,fileName,num,uce,inuce,window,nucLine,methFlank)
	if 'dendrogram' in graphs:
		GraphDendrogramLibrary.main(dfWindow,ranWindow,names,fileName,num,uce,inuce,window,nucLine,methFlank)
	if 'kmean' in graphs:
		GraphKMeansLibrary.main(dfWindow,ranWindow,names,fileName,num,uce,inuce,window,nucLine,methFlank)

# make a table for element and random region
def plotTable(TableData,Title,ranTableData,ranTitle,fileName):
	GraphTableLibrary.main(TableData,Title,ranTableData,ranTitle,fileName)

# collect the counts for how many of each boundary direction for each type
def collectCounts(rangeFeatures):
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
	
	# Integer parameters
	num = args.total
	uce = args.element
	inuce = args.inset
	window = args.window
	binDir = args.bin
	methCovThresh = args.thresholdcoverage
	methPerThresh = args.thresholdpercentage
	methFlank = args.methylationflank

	# Element, random regions and methylation files
	eFiles = [line.strip() for line in args.efile]
	mFiles = [line.strip() for line in args.mfile]
	rFiles = [line.strip() for line in args.rfile]

	# Genome files from UCSC
	sizeGenome = args.genome
	faGenome = args.fasta
	Overlapregions = args.overlapingelements

	# Lists with the types and directions to use
	typeList = args.elementype
	dirList = args.elementdirection
	nucLine = args.nucleotideline

	# Reverse complement argument
	revCom = args.reversecomplement
	overlapInset = args.elementalign

	# Which plots to run
	graphs = args.plots

	# A string to add to the out file name in case you want to set up runs and let be
	stringName = args.stringname

	# for each element file provided
	for fileName in eFiles:
		# get the region info to work with
		rangeFeatures = ElementLibrary.main(num,uce,inuce,window,binDir,fileName,sizeGenome,faGenome)
		# separate by direction
		directionFeatures = DirectionLibrary.main(rangeFeatures,fileName,binDir)
		
		# for each random file provided
		for randomFile in rFiles:
			# Collect out filename labels 
			paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(uce,inuce,num,binDir,window,fileName,randomFile,stringName)
			
			# get the region info to work with
			randomFeatures = ElementLibrary.main(num,uce,inuce,window,binDir,randomFile,sizeGenome,faGenome)
			# separate by direction
			randirFeatures = DirectionLibrary.main(randomFeatures,randomFile,binDir)
			
			# Make table for the count of each direction for each type element
			elementGroup = collectCounts(rangeFeatures)
			randomGroup = collectCounts(randirFeatures)
			plotTable(elementGroup,'Elements',randomGroup,'Random Regions','Direction_Count_{0}'.format(paramlabels))
			
			# All elements
			if 'all' in typeList:
				typeList.remove('all')
				pdMeth,allWindow, allNames = TypeLibrary.main(rangeFeatures,'combineString',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
				rnMeth,ranWindow, ranNames = TypeLibrary.main(randomFeatures,'combineString',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
				plotGraphs(pdMeth,rnMeth,allWindow,allNames,ranWindow,'all_{0}'.format(paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)
				if revCom:
					revMeth,revWindow,revNames = RevCompLibrary.main(directionFeatures,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
					ranrevMeth,ranrevWindow,ranrevNames = RevCompLibrary.main(randirFeatures,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
					plotGraphs(revMeth,ranrevMeth,revWindow,revNames,ranrevWindow,'revComp_all_{0}'.format(paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)

# 			# By Type
			for type in typeList:
				typeBool,typeMeth,typeWindow,typeNames = groupSeparate(type,directionFeatures,'type',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs,'combineString')
				rantypeBool,rantypeMeth,rantypeWindow,rantypeNames = groupSeparate(type,randirFeatures,'type',randomFile,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs,'combineString')
				plotGraphs(typeMeth,rantypeMeth,typeWindow,typeNames,rantypeWindow,'{0}_{1}'.format(type,paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)
				if revCom:
					typercMeth,typercWindow,typercNames = RevCompLibrary.main(typeBool,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
					rantypercMeth,rantypercWindow,rantypercNames = RevCompLibrary.main(rantypeBool,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
					plotGraphs(typercMeth,rantypercMeth,typercWindow,typercNames,rantypercWindow,'revComp_{0}_{1}'.format(type,paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)

			# By Direction
			for dir in dirList:
				dirBool,dirMeth,dirWindow,dirNames = groupSeparate(dir,directionFeatures,'compareBoundaries',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs,'combineString')
				randirBool,randirMeth,randirWindow,randirNames = groupSeparate(dir,randirFeatures,'compareBoundaries',randomFile,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs,'combineString')
				plotGraphs(dirMeth,randirMeth,dirWindow,dirNames,randirWindow,'all_{0}_{1}'.format(dir,paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)
				for type in typeList:
					typeBool,typeMeth,typeWindow,typeNames = groupSeparate(type,directionFeatures,'type',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs,'combineString')
					rantypeBool,rantypeMeth,rantypeWindow,rantypeNames = groupSeparate(type,randirFeatures,'type',randomFile,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs,'combineString')
					plotGraphs(typeMeth,rantypeMeth,typeWindow,typeNames,rantypeWindow,'{0}_{1}_{2}'.format(type,dir,paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)

			# Re-align by exon/intron crossover
			if overlapInset:
				startcrossboundary,endcrossboundary,completeelement,interiorelement = OverlapLibrary.main(rangeFeatures,Overlapregions,num,uce,inuce,faGenome,binDir,revCom,fileName,mFiles,window,methCovThresh,methPerThresh,nucLine,graphs)
				crossMeth,crossWindow,crossNames  = RevCompOverlapLibrary.main(startcrossboundary,endcrossboundary,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
				completeelement.name = 'WithinExon'
				interiorelement.name = 'ContainsExon'
				ranstartcrossboundary,ranendcrossboundary,rancompleteelement,raninteriorelement = OverlapLibrary.main(randomFeatures,Overlapregions,num,uce,inuce,faGenome,binDir,revCom,fileName,mFiles,window,methCovThresh,methPerThresh,nucLine,graphs)
				rancrossMeth,rancrossWindow,rancrossNames  = RevCompOverlapLibrary.main(ranstartcrossboundary,ranendcrossboundary,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
				plotGraphs(crossMeth,rancrossMeth,crossWindow,crossNames,rancrossWindow,'align_Crossboundary_{0}'.format(paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)
				for element,random in zip([completeelement,interiorelement],[rancompleteelement,raninteriorelement]):
					alignMeth,alignWindow,alignNames = TypeLibrary.main(element,'combineString',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
					ranalignMeth,ranalignWindow,ranalignNames = TypeLibrary.main(random,'combineString',fileName,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
					plotGraphs(alignMeth,ranalignMeth,alignWindow,alignNames,ranalignWindow,'align_{0}_{1}'.format(element.name,paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)
# 					if revCom: # maybe use group separate?
# 						aligndirFeatures = DirectionLibrary.main(element,fileName,binDir)
# 						alignranFeatures = DirectionLibrary.main(random,fileName,binDir)
# 						alignrcMeth,alignrcWindow,alignrcNames = RevCompLibrary.main(aligndirFeatures,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
# 						ranalignrcMeth,ranalignrcWindow,ranalignrcNames = RevCompLibrary.main(alignranFeatures,binDir,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
# 						plotGraphs(alignrcMeth,ranalignrcMeth,alignrcWindow,alignrcNames,ranalignrcWindow,'align_{0}_{1}'.format(element.name,paramlabels),num,uce,inuce,window,graphs,nucLine,methFlank)





if __name__ == "__main__":
	main()