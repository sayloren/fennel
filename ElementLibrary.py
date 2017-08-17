"""
Script to process the elements you what to analyze

Wren Saylor
July 5 2017
"""

import argparse
from collections import defaultdict
import pandas as pd
import pybedtools as pbt
import GlobalVariables

# read in files
def eachFileProcess(fileName):
	btFeatures = pbt.BedTool(fileName)
	return btFeatures

# get bt features
def getFeatures(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# get the correct range for fang evaluation
def getRange(btFeatures,fileName):#btFeatures,fileName,num,uce,inuce
	flankSize = (GlobalVariables.num - GlobalVariables.uce)/2
	inregion = GlobalVariables.uce-(GlobalVariables.inuce*2)
	midFeatures = pd.read_table(btFeatures.fn, header=None)
	midFeatures['middle'] = midFeatures.loc[:,1:2].mean(axis=1).astype(int)
	midFeatures['start'] = midFeatures.loc[:,1]
	midFeatures['end'] = midFeatures.loc[:,2] 
	midFeatures['sCenter'] = midFeatures['middle'] - (inregion/2)
	midFeatures['eCenter'] = midFeatures['middle'] + (inregion/2)
	midFeatures['sEdge'] = midFeatures.loc[:,1] + GlobalVariables.inuce
	midFeatures['eEdge'] = midFeatures.loc[:,2] - GlobalVariables.inuce
	midFeatures['sBoundary'] = midFeatures.loc[:,1] - flankSize
	midFeatures['eBoundary'] = midFeatures.loc[:,2] + flankSize
	midFeatures['type'] = midFeatures.loc[:,4]
	midFeatures['id'] = midFeatures.loc[:,3]
	midFeatures['chr'] = midFeatures.loc[:,0]
	midFeatures['size'] = midFeatures.loc[:,2].astype(int)-midFeatures.loc[:,1].astype(int)
	rangeFeatures = midFeatures[['type','id','size','chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary']]
	print 'Getting the coordinates for the area to examine, {0} from the middle of the element, {1} inset from the edges'.format((midFeatures['eCenter']-midFeatures['sCenter']).mean(),(midFeatures['sEdge']-midFeatures['start']).mean())
	return rangeFeatures

# get the strings for sliding window regions
def btRange(rangeFeatures):#rangeFeatures,faGenome
	rangeFeatures['sBoundarySeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','sBoundary','start']].values.tolist()))
	rangeFeatures['sEdgeSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','start','sEdge']].values.tolist()))
	rangeFeatures['MiddleSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','sCenter','eCenter']].values.tolist()))
	rangeFeatures['eEdgeSeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','eEdge','end',]].values.tolist()))
	rangeFeatures['eBoundarySeq'] = simpleFasta(getFeatures(rangeFeatures[['chr','end','eBoundary']].values.tolist()))
	rangeFeatures['feature'] = simpleFasta(getFeatures(rangeFeatures[['chr','start','end']].values.tolist()))
	rangeFeatures['combineString'] = rangeFeatures['sBoundarySeq'].astype(str) + rangeFeatures['sEdgeSeq'].astype(str) + rangeFeatures['MiddleSeq'].astype(str) + rangeFeatures['eEdgeSeq'].astype(str) + rangeFeatures['eBoundarySeq'].astype(str)
	rangeFeatures['combineString'] = rangeFeatures['combineString'].str.upper()
	rangeFeatures['feature'] = rangeFeatures['feature'].str.upper()
	rangeFeatures['reverseComplement'] = rangeFeatures.apply(lambda row: reverseComplement(row['combineString']),axis=1)
	return rangeFeatures

# Get the percentage AT in the element
def perElementAT(region,fileName):
	collectAT = []
	for r in region:
		collectAT.append(eval('100*float(r.count("A") + r.count("a") + r.count("T") + r.count("t"))/len(r)'))
	pdAT = pd.DataFrame(collectAT)
	print 'Mean AT content for all {0} elements in {1} is {2} Percent'.format(len(region.index),fileName,pdAT.mean())

# get the reverse complement
def reverseComplement(sequence):
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	return "".join([seqDict[base] for base in reversed(sequence)])

# save file from bedtool
def saveBedTool(btObject,strFilename):
	btObject.saveas(strFilename)

# convert bedtool to panda
# If there is nothing in the btobject, will it read the data from the previous itteration!?
def bedtoolToPanda(btobject):
	saveBedTool(btobject,'temp.bed')
	pdObject = pd.read_table(btobject.fn, header=None)
	return pdObject

# convert panda to bedtool
def pandaToBedtool(panda):
	arArFeatures = panda.values.tolist()
	btoutFeatures = getFeatures(arArFeatures)
	return btoutFeatures

# get fasta strings for each desired region
def getFasta(btFeatures,fileName):
	saveAll = 'Seq_result_for_all_{0}.txt'.format(fileName)
	seqAll = btFeatures.sequence(fi=GlobalVariables.faGenome,fo=saveAll)
	saveExonic = 'Seq_results_for_exonic_{0}.txt'.format(fileName)
	seqExonic = btFeatures.sequence(fi=GlobalVariables.faGenome,fo=saveExonic).filter(lambda x: x[name] == 'exonic')
	saveIntronic = 'Seq_results_for_intronic_{0}.txt'.format(fileName)
	seqIntronic = btFeatures.sequence(fi=GlobalVariables.faGenome,fo=saveIntronic).filter(lambda x: x[name] == 'intronic')
	saveIntergenic = 'Seq_results_for_intergenic_{0}.txt'.format(fileName)
	seqIntergenic = btFeatures.sequence(fi=vfaGenome,fo=saveIntergenic).filter(lambda x: x[name] == 'intergenic')
	return saveAll, saveExonic, saveIntronic, saveIntergenic

# used in btRange to extract just the fasta strings
def simpleFasta(inFeature):
	seqFeature = inFeature.sequence(fi=GlobalVariables.faGenome)
	outFeature = pd.read_table(seqFeature.seqfn)
	outSequence = outFeature[::2]
	outSequence = outSequence.reset_index(drop=True)
	return outSequence

def main(fileName):
	print 'Running ElementLibrary'
	btFeatures = eachFileProcess(fileName)
	subsetFeatures = getRange(btFeatures,fileName)
	rangeFeatures = btRange(subsetFeatures)
	perElementAT(rangeFeatures['feature'],fileName)
	return rangeFeatures

if __name__ == "__main__":
	main()