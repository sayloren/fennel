"""
Script to align exons

Wren Saylor
August 4 2017
"""

import argparse
from collections import defaultdict
import pandas as pd
import pybedtools as pbt
from ElementLibrary import getRange
from ElementLibrary import btRange
from ElementLibrary import eachFileProcess
from ElementLibrary import bedtoolToPanda
import DirectionLibrary
from RevCompLibrary import slideDirection

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

def rangebyClass(df,num,uce,inuce):
	flankSize = (num - uce)/2
	inregion = uce-(inuce*2)
# 	start = df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)].count() # how many elements overlap at the upstream boundary
# 	end = df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)].count() # how many elements overlap at the downstream boundary
# 	interior = df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)].count() # how many elements contain an exon
# 	multi = df[df.duplicated(subset='id',keep='first')].count() # number of elements that fit into two categories
	if df[(df['startdifference'] >= 0) & (df['enddifference']  >= 0)] # elements are completely inside the exon
		midFeatures['middle'] = midFeatures[['start','end']].mean(axis=1).astype(int)
		midFeatures['sCenter'] = midFeatures['middle'] - (inregion/2)
		midFeatures['eCenter'] = midFeatures['middle'] + (inregion/2)
		midFeatures['sEdge'] = midFeatures[['start']] + inuce
		midFeatures['eEdge'] = midFeatures[['end']] - inuce
		midFeatures['sBoundary'] = midFeatures[['start']] - flankSize
		midFeatures['eBoundary'] = midFeatures[['end']] + flankSize
	else:
		out = None




# 	midFeatures['middle'] = midFeatures[['start','end']].mean(axis=1).astype(int)
# 	midFeatures['sCenter'] = midFeatures['middle'] - (inregion/2)
# 	midFeatures['eCenter'] = midFeatures['middle'] + (inregion/2)
# 	midFeatures['sEdge'] = midFeatures[['start']] + inuce
# 	midFeatures['eEdge'] = midFeatures[['end']] - inuce
# 	midFeatures['sBoundary'] = midFeatures[['start']] - flankSize
# 	midFeatures['eBoundary'] = midFeatures[['end']] + flankSize
# 	rangeFeatures = midFeatures[['id','chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary']]
# 	return rangeFeatures
# 
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
	rangebyClass(exonicFeatures,num,uce,inuce)

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