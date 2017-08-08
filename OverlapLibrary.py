"""
Script to align exons

Wren Saylor
August 4 2017
"""

import argparse
from collections import defaultdict
import pandas as pd
import pybedtools as pbt
from ElementLibrary import reverseComplement
from ElementLibrary import eachFileProcess
from ElementLibrary import bedtoolToPanda
from ElementLibrary import simpleFasta
from ElementLibrary import getFeatures
# import DirectionLibrary
# from RevCompLibrary import slideDirection

# Intersect the exons with the uces
def exonIntersect(df,exon,num,uce,inuce,faGenome):
	ExonBoundary = exon.intersect(df[['chr','start','end','id']].values.tolist(),wo=True)
	pdExonBoundary = bedtoolToPanda(ExonBoundary)
	pdExonBoundary.columns = ['ochr','ostart','oend','chr','start','end','id','overlap']#,'exondir'
	pdExonBoundary['elementsize'] = pdExonBoundary['end'] - pdExonBoundary['start']
	pdExonBoundary['startdifference'] = pdExonBoundary['start'] - pdExonBoundary['ostart']# if this is a negative number, the element starts before the exon
	pdExonBoundary['enddifference'] = pdExonBoundary['oend'] - pdExonBoundary['end']# if this is a negative number, the element ends after the exon
	overlapTable = makeOverlaptable(pdExonBoundary)
	print overlapTable
	return pdExonBoundary, overlapTable

def makeOverlaptable(df):
	complete = df[(df['startdifference'] >= 0) & (df['enddifference']  >= 0)].count() # how many elements are completely inside the exon
	start = df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)].count() # how many elements overlap at the upstream boundary
	end = df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)].count() # how many elements overlap at the downstream boundary
	interior = df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)].count() # how many elements contain an exon
	multi = df[df.duplicated(subset='id',keep='first')].count() # number of elements that fit into two categories
	total = df.count()-multi # get the total number of elements intersecting exons - minus those that overlap two exons
	frames = [complete,start,end,interior,multi,total]
	overlapTable = pd.concat(frames,axis=1)
	overlapTable.columns = ['complete','upstream','downstream','interior','multiple','total']
	overlapTable = overlapTable.head(1)
	overlapTable.index = ['overlaps']
	return overlapTable

def rangebyClass(df,num,uce,inuce,faGenome): #how will rc sorting work? maybe just for complete and interior
	centerElement = num/2
	flankSize = (num - uce)/2
	inregion = uce-(inuce*2)

	if len(df[(df['startdifference'] >= 0) & (df['enddifference'] >= 0)]) != 0: # elements are completely inside the exon
		features = df[(df['startdifference'] >= 0) & (df['enddifference'] >= 0)]
		features['middle'] = features[['start','end']].mean(axis=1).astype(int)
		features['sCenter'] = features['middle'] - (inregion/2)
		features['eCenter'] = features['middle'] + (inregion/2)
		features['sEdge'] = features[['start']] + inuce
		features['eEdge'] = features[['end']] - inuce
		features['sBoundary'] = features[['start']] - flankSize
		features['eBoundary'] = features[['end']] + flankSize
		features['sBoundarySeq'] = simpleFasta(getFeatures(features[['chr','sBoundary','start']].values.tolist()),faGenome)
		features['sEdgeSeq'] = simpleFasta(getFeatures(features[['chr','start','sEdge']].values.tolist()),faGenome)
		features['MiddleSeq'] = simpleFasta(getFeatures(features[['chr','sCenter','eCenter']].values.tolist()),faGenome)
		features['eEdgeSeq'] = simpleFasta(getFeatures(features[['chr','eEdge','end',]].values.tolist()),faGenome)
		features['eBoundarySeq'] = simpleFasta(getFeatures(features[['chr','end','eBoundary']].values.tolist()),faGenome)
		features['feature'] = simpleFasta(getFeatures(features[['chr','start','end']].values.tolist()),faGenome)
		features['combineString'] = features['sBoundarySeq'].astype(str) + features['sEdgeSeq'].astype(str) + features['MiddleSeq'].astype(str) + features['eEdgeSeq'].astype(str) + features['eBoundarySeq'].astype(str)
		features['combineString'] = features['combineString'].str.upper()
		features['reverseComplement'] = features.apply(lambda row: reverseComplement(row['combineString']),axis=1) # may use once combine features?
		completeelement = features[['chr','sBoundary','eBoundary','combineString','reverseComplement','id']]
		completeelement.columns = ['chr','start','end','string','revstring','id']
	else:
		completeelement = None
		
	if len(df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)]) != 0: # elements contain an exon
		features = df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)]
		interiorsize = df['overlap'].min() # smallest exon, to be able to compare the exon boundaries crossover
		df['difference'] = df['elementsize'] - df['overlap']
		elementsize = df['difference'].min()
		features['sCenter'] = features['ostart'] + (interiorsize/2)
		features['eCenter'] = features['oend'] - (interiorsize/2)
# 		features['usEdge'] = features[['start']] # upstream
# 		features['ueEdge'] = features[['start']] + (elementsize/2) # upstream
# 		features['dsEdge'] = features[['end']] # downstream
# 		features['deEdge'] = features[['end']]  - (elementsize/2) # downstream
		features['sBoundary'] = features[['start']] - flankSize
		features['eBoundary'] = features[['end']] + flankSize
		
		
		features['sBoundarySeq'] = simpleFasta(getFeatures(features[['chr','sBoundary','start']].values.tolist()),faGenome)
# 		features['usEdgeSeq'] = simpleFasta(getFeatures(features[['chr','start','sEdge']].values.tolist()),faGenome)
# 		features['ueEdgeSeq'] = simpleFasta(getFeatures(features[['chr','start','sEdge']].values.tolist()),faGenome)
		features['MiddleSeq'] = simpleFasta(getFeatures(features[['chr','sCenter','eCenter']].values.tolist()),faGenome)
# 		features['dsEdgeSeq'] = simpleFasta(getFeatures(features[['chr','eEdge','end',]].values.tolist()),faGenome)
# 		features['deEdgeSeq'] = simpleFasta(getFeatures(features[['chr','eEdge','end',]].values.tolist()),faGenome)
		features['eBoundarySeq'] = simpleFasta(getFeatures(features[['chr','end','eBoundary']].values.tolist()),faGenome)
		features['feature'] = simpleFasta(getFeatures(features[['chr','start','end']].values.tolist()),faGenome)
		features['combineString'] = features['sBoundarySeq'].astype(str) + features['sEdgeSeq'].astype(str) + features['MiddleSeq'].astype(str) + features['eEdgeSeq'].astype(str) + features['eBoundarySeq'].astype(str)
		features['combineString'] = features['combineString'].str.upper()
		features['reverseComplement'] = features.apply(lambda row: reverseComplement(row['combineString']),axis=1) # may use once combine features?

		
		interiorelement = features[[]]
		interiorelement.columns = []
	else:
		interiorelement = None

# Exonic cross boundary shift (still need to separate by intergenic/intronic)
	if len(df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)]) != 0: # elements overlap at the upstream boundary
		features = df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)]
		features['startCoord'] = features[['oend']] - centerElement
		features['crossBoundary'] = features[['oend']]
		features['endCoord'] = features[['oend']] + centerElement
		features['nucString'] = simpleFasta(getFeatures(features[['chr','startCoord','endCoord']].values.tolist()),faGenome)
		features['nucString'] = features['nucString'].str.upper()
		startout = features[['chr','startCoord','endCoord','nucString','id']]
		startout.columns = ['chr','start','end','string','id']
	else:
		startout = None

	if len(df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)]) != 0: # elements overlap at the downstream boundary
		features = df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)]
		features['startCoord'] = features[['ostart']] - centerElement
		features['crossBoundary'] = features[['ostart']]
		features['endCoord'] = features[['ostart']] + centerElement
		features['nucString'] = simpleFasta(getFeatures(features[['chr','startCoord','endCoord']].values.tolist()),faGenome)
		features['nucString'] = features['nucString'].str.upper()
		features['reverseComplement'] = features.apply(lambda row: reverseComplement(row['nucString']),axis=1)
		endout = features[['chr','startCoord','endCoord','reverseComplement','id']]
		endout.columns = ['chr','start','end','string','id']
	else:
		endout = None
	boundaryframes = [startout,endout]
	crossboundary = pd.concat(boundaryframes)

	
	return crossboundary, completeelement, interiorelement
	
	
	
# def mostlydirLine(directionFeatures,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs):
# 	negStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '-')])
# 	posStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '+')])
# 	negCoord = mostlyGetRange(negStr,num,uce,inuce)
# 	posCoord = mostlyGetRange(posStr,num,uce,inuce)
# 	negFeat = btRange(negCoord,faGenome)
# 	posFeat = btRange(posCoord,faGenome)
# 	compWindow, compNames = slideDirection(negFeat,posFeat,num,uce,inuce,window,nucLine)
# 	if any(x in graphs for x in ['methylation','cluster']):
# 		groupMeth = methDirection(negFeat,posFeat,mFiles,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
# 	else: 
# 		groupMeth = None
# 	return groupMeth,compWindow,compNames

def main(rangeFeatures,exonicInset,num,uce,inuce,faGenome,binDir,revCom,fileName,mFiles,window,methCovThresh,methPerThresh,nucLine,graphs):
	exonFeature = eachFileProcess(exonicInset)
	exonicFeatures, overlapTable = exonIntersect(rangeFeatures,exonFeature,num,uce,inuce,faGenome)
	rangebyClass(exonicFeatures,num,uce,inuce,faGenome)

# 	subsetFeatures = mostlyGetRange(exonicFeatures,num,uce,inuce)
# 	rangeFeatures = btRange(subsetFeatures,faGenome)
# 	directionFeatures = DirectionLibrary.main(rangeFeatures,fileName,binDir)
# 	if revCom:
# 		groupMeth,compWindow,compNames = mostlydirLine(directionFeatures,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
# 	else: 
# 		rcdirectionFeatures = None
# 	return directionFeatures, rcdirectionFeatures

if __name__ == "__main__":
	main()