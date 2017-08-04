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
def exonIntersect(rangeFeatures,exonFeature,num,uce,inuce,faGenome):
	ExonBoundary = exonFeature.intersect(rangeFeatures[['chr','start','end','id']].values.tolist(),wb=True,wa=True)
	pdExonBoundary = bedtoolToPanda(ExonBoundary)
	pdExonBoundary.columns = ['echr','estart','estop','chr','start','end','id']#,'exondir'
	subExonBoundary = pdExonBoundary[['echr','estart','end','id']]#,'exondir'
	subExonBoundary.columns = ['chr','start','end','id']#,'exondir'
	return subExonBoundary

def mostlyGetRange(midFeatures,num,uce,inuce):
	flankSize = (num - uce)/2
	inregion = uce-(inuce*2)
	midFeatures['middle'] = midFeatures[['start','end']].mean(axis=1).astype(int)
	midFeatures['sCenter'] = midFeatures['middle'] - (inregion/2)
	midFeatures['eCenter'] = midFeatures['middle'] + (inregion/2)
	midFeatures['sEdge'] = midFeatures[['start']] + inuce
	midFeatures['eEdge'] = midFeatures[['end']] - inuce
	midFeatures['sBoundary'] = midFeatures[['start']] - flankSize
	midFeatures['eBoundary'] = midFeatures[['end']] + flankSize
	rangeFeatures = midFeatures[['id','chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary']]
	return rangeFeatures

def mostlydirLine(directionFeatures,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs):
	negStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '-')])
	posStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '+')])
	negCoord = mostlyGetRange(negStr,num,uce,inuce)
	posCoord = mostlyGetRange(posStr,num,uce,inuce)
	negFeat = btRange(negCoord,faGenome)
	posFeat = btRange(posCoord,faGenome)
	compWindow, compNames = slideDirection(negFeat,posFeat,num,uce,inuce,window,nucLine)
	if any(x in graphs for x in ['methylation','cluster']):
		groupMeth = methDirection(negFeat,posFeat,mFiles,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
	else: 
		groupMeth = None
	return groupMeth,compWindow,compNames

def main(rangeFeatures,exonicInset,num,uce,inuce,faGenome,binDir,revCom,fileName,mFiles,window,methCovThresh,methPerThresh,nucLine,graphs):
	exonFeature = eachFileProcess(exonicInset)
	exonicFeatures = exonIntersect(rangeFeatures,exonFeature,num,uce,inuce,faGenome)
	subsetFeatures = mostlyGetRange(exonicFeatures,num,uce,inuce)
	rangeFeatures = btRange(subsetFeatures,faGenome)
	directionFeatures = DirectionLibrary.main(rangeFeatures,fileName,binDir)
	if revCom:
		groupMeth,compWindow,compNames = mostlydirLine(directionFeatures,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
	else: 
		rcdirectionFeatures = None
	return directionFeatures, rcdirectionFeatures

if __name__ == "__main__":
	main()