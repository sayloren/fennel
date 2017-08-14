"""
Script to separate by directionality

Wren Saylor
July 5 2017
"""
import argparse
import numpy as np
import pandas as pd

def calculateAT(element,size):
	start = element[:size]
	end = element[-size:]
	perSize = []
	perSize.append(eval('100*float(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
	perSize.append(eval('100*float(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
	return perSize

# out put directionality, as inferred by comparing first and last n base pairs from input parameters
def compareN(element,size):
	perSize = calculateAT(element,size)
	# give + - = depending on which side has larger AT content
	if perSize[0] > perSize[1]: outList = '+'
	if perSize[1] > perSize[0]: outList = '-'
	if perSize[1] == perSize[0]: outList = '='
	return outList

# with the results from compareN per each element, evaluate directionality into new column
def evalN(rangeFeatures,fileName,binDir):
	rangeFeatures['compareBoundaries'] = rangeFeatures.apply(lambda row: (compareN(row['feature'],binDir)),axis=1)
	compareEnds = pd.DataFrame(rangeFeatures[['chr','start','end','compareBoundaries']])
	# put bin size calibration here
	print 'Sorting the element boundaries by bin size {0}'.format(binDir)
	
	rangeBins = collectEmperical(rangeFeatures,binDir)
	# for after rc sorting too?
	return rangeFeatures,rangeBins

def collectEmperical(rangeFeatures,binDir):
	# Perform min,max collection
	rangeAT = rangeFeatures.apply(lambda row: (empericalATspread(row['feature'],binDir)),axis=1)
	pdAT = pd.DataFrame(rangeAT.values.tolist())
	
	# Separate Start from End
	splitAT = pd.concat(dict([(row[0],row[1].apply(lambda y: pd.Series(y))) for row in pdAT.iterrows()]),axis=1)
	
	outcollect = []
	for column in splitAT:
		outcollect.append(splitAT[column[0]])
	outcat = pd.concat(outcollect,axis=1)
	outcat /= 100 # convert to decimal

	# Separate data frames by Start v End and Min v Max
	pdATCollectStartMin = outcat[[outcat.columns[0]]].min(axis=1)
	pdATCollectEndMin = outcat[[outcat.columns[1]]].min(axis=1)
	pdATCollectStartMax = outcat[[outcat.columns[0]]].max(axis=1)
	pdATCollectEndMax = outcat[[outcat.columns[1]]].max(axis=1)

	Min = pdATCollectStartMin * pdATCollectEndMin
	Max = pdATCollectStartMax * pdATCollectEndMax

	pdBins = pd.concat([Min,Max],axis=1)
	pdBins.columns=['Min','Max']
	return pdBins

def empericalATspread(element,size):
	totalSteps=[]
	for i in np.arange(1,size*2):
		pairStep = calculateAT(element,i)
		totalSteps.append(pairStep)
	return totalSteps

def main(rangeFeatures,fileName,binDir):
	directionFeatures,directionBins = evalN(rangeFeatures,fileName,binDir)
	return directionFeatures,directionBins

if __name__ == "__main__":
	main()