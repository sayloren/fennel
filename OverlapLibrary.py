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

def savePanda(pdData, strFilename):
	pdData.to_csv(strFilename, sep='\t')

# Intersect the exons with the uces
def exonIntersect(df,exon,num,uce,inuce,faGenome):
	ExonBoundary = exon.intersect(df[['chr','start','end','id']].values.tolist(),wo=True)
	pdExonBoundary = bedtoolToPanda(ExonBoundary)
	pdExonBoundary.columns = ['ochr','ostart','oend','chr','start','end','id','overlap']#,'exondir'
	pdExonBoundary['elementsize'] = pdExonBoundary['end'] - pdExonBoundary['start']
	pdExonBoundary['startdifference'] = pdExonBoundary['start'] - pdExonBoundary['ostart']# if this is a negative number, the element starts before the exon
	pdExonBoundary['enddifference'] = pdExonBoundary['oend'] - pdExonBoundary['end']# if this is a negative number, the element ends after the exon
	overlapTable = makeOverlaptable(pdExonBoundary)
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

def rangebyClass(df,num,uce,inuce,faGenome):
	centerElement = num/2
	flankSize = (num - uce)/2
	inregion = uce-(inuce*2)
	# Exonic regions that are completely within an exon
	if len(df[(df['startdifference'] >= 0) & (df['enddifference'] >= 0)]) != 0: # elements are completely inside the exon
		features = df[(df['startdifference'] >= 0) & (df['enddifference'] >= 0)]
		features.reset_index(drop=True,inplace=True)
		features['middle'] = features[['start','end']].mean(axis=1).astype(int)
		features['sCenter'] = features.loc[:,'middle'] - (inregion/2)
		features['eCenter'] = features.loc[:,'middle'] + (inregion/2)
		features['sEdge'] = features.loc[:,'start'] + inuce
		features['eEdge'] = features.loc[:,'end'] - inuce
		features['sBoundary'] = features.loc[:,'start'] - flankSize
		features['eBoundary'] = features.loc[:,'end'] + flankSize
		features['sBoundarySeq'] = simpleFasta(getFeatures(features[['chr','sBoundary','sEdge']].values.tolist()),faGenome)
		features['MiddleSeq'] = simpleFasta(getFeatures(features[['chr','sCenter','eCenter']].values.tolist()),faGenome)
		features['eBoundarySeq'] = simpleFasta(getFeatures(features[['chr','eEdge','eBoundary']].values.tolist()),faGenome)
		features['feature'] = simpleFasta(getFeatures(features[['chr','start','end']].values.tolist()),faGenome)
		features['combineString'] = features['sBoundarySeq'].astype(str) +features['MiddleSeq'].astype(str) + features['eBoundarySeq'].astype(str)
		features['combineString'] = features['combineString'].str.upper()
		features['size'] = features['combineString'].str.len() # check length of string
		features['reverseComplement'] = features.apply(lambda row: reverseComplement(row['combineString']),axis=1) # may use once combine features?
		completeelement = features[['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','id']]
		completeelement.columns = ['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','id']
	else:
		completeelement = None
	# Exonic regions that contain an exon - was going to align by interior exon, buuuut, its tricky
	if len(df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)]) != 0: # elements contain an exon
		features = df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)]
		interiorsize = df['overlap'].min() # smallest exon, to be able to compare the exon boundaries crossover
		df['difference'] = df['elementsize'] - df['overlap']
		elementsize = df['difference'].min()
		features = df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)]
		features.reset_index(drop=True,inplace=True)
		features['middle'] = features[['start','end']].mean(axis=1).astype(int)
		features['sCenter'] = features['middle'] - (inregion/2)
		features['eCenter'] = features['middle'] + (inregion/2)
		features['sEdge'] = features['start'] + inuce
		features['eEdge'] = features['end'] - inuce
		features['sBoundary'] = features['start'] - flankSize
		features['eBoundary'] = features['end'] + flankSize
		features['sBoundarySeq'] = simpleFasta(getFeatures(features[['chr','sBoundary','start']].values.tolist()),faGenome)
		features['sEdgeSeq'] = simpleFasta(getFeatures(features[['chr','start','sEdge']].values.tolist()),faGenome)
		features['MiddleSeq'] = simpleFasta(getFeatures(features[['chr','sCenter','eCenter']].values.tolist()),faGenome)
		features['eEdgeSeq'] = simpleFasta(getFeatures(features[['chr','eEdge','end',]].values.tolist()),faGenome)
		features['eBoundarySeq'] = simpleFasta(getFeatures(features[['chr','end','eBoundary']].values.tolist()),faGenome)
		features['feature'] = simpleFasta(getFeatures(features[['chr','start','end']].values.tolist()),faGenome)
		features['combineString'] = features['sBoundarySeq'].astype(str) + features['sEdgeSeq'].astype(str) + features['MiddleSeq'].astype(str) + features['eEdgeSeq'].astype(str) + features['eBoundarySeq'].astype(str)
		features['combineString'] = features['combineString'].str.upper()
		features['size'] = features['combineString'].str.len() # check length of string
		features['reverseComplement'] = features.apply(lambda row: reverseComplement(row['combineString']),axis=1) # may use once combine features?
		interiorelement = features[['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','id']]
		interiorelement.columns = ['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','id']
	else:
		interiorelement = None
	# Exonic cross boundary shift (still need to separate by intergenic/intronic)
# 	if len(df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)]) != 0: # elements overlap at the upstream boundary
# 		features = df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)]
# 		features.reset_index(drop=True,inplace=True)
# 		features['sEdge'] = features['oend']
# 		features['sCenter'] = features['oend'] + (inregion/2)
# 		features['eCenter'] = features['oend'] + uce - (inregion/2)
# 		features['eEdge'] = features['oend'] + uce
# 		features['sBoundary'] = features['oend'] - flankSize
# 		features['eBoundary'] = features['eEdge'] + flankSize
# 		features['combineString'] = simpleFasta(getFeatures(features[['chr','sBoundary','eBoundary']].values.tolist()),faGenome)
# 		features['combineString'] = features['combineString'].str.upper()
# 		features['reverseComplement'] = features.apply(lambda row: reverseComplement(row['combineString']),axis=1) # may use once combine features?
# 
# 		features['size'] = features['combineString'].str.len() # check length of string
# 		startout = features[['chr','startCoord','endCoord','combineString','id']]
# 		startout.columns = ['chr','start','end','combineString','id']
# 	else:
# 		startout = None
# 	if len(df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)]) != 0: # elements overlap at the downstream boundary
# 		features = df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)]
# 		features.reset_index(drop=True,inplace=True)
# 		features['sEdge'] = features['ostart']
# 		features['sCenter'] = features['ostart'] + (inregion/2)
# 		features['eCenter'] = features['ostart'] + uce - (inregion/2)
# 		features['eEdge'] = features['ostart'] + uce
# 		features['sBoundary'] = features['ostart'] - flankSize
# 		features['eBoundary'] = features['eEdge'] + flankSize
# 		features['combineString'] = simpleFasta(getFeatures(features[['chr','sBoundary','eBoundary']].values.tolist()),faGenome)
# 		features['combineString'] = features['combineString'].str.upper()
# 		features['reverseComplement'] = features.apply(lambda row: reverseComplement(row['combineString']),axis=1) # may use once combine features?

		features['size'] = features['combineString'].str.len() # check length of string
		endout = features[['chr','startCoord','endCoord','reverseComplement','id']]
		endout.columns = ['chr','start','end','combineString','id']
	else:
		endout = None
	boundaryframes = [startout,endout]
	crossboundary = pd.concat(boundaryframes)
	crossboundary.reset_index(drop=True,inplace=True)
	return crossboundary,completeelement,interiorelement

def main(rangeFeatures,exonicInset,num,uce,inuce,faGenome,binDir,revCom,fileName,mFiles,window,methCovThresh,methPerThresh,nucLine,graphs):
	exonFeature = eachFileProcess(exonicInset)
	exonicFeatures, overlapTable = exonIntersect(rangeFeatures,exonFeature,num,uce,inuce,faGenome)
	crossboundary,completeelement,interiorelement = rangebyClass(exonicFeatures,num,uce,inuce,faGenome)
	return crossboundary,completeelement,interiorelement

if __name__ == "__main__":
	main()