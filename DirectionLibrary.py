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
	perSize.append(eval('100*int(start.count("A") + start.count("a") + start.count("T") + start.count("t"))/len(start)'))
	perSize.append(eval('100*int(end.count("A") + end.count("a") + end.count("T") + end.count("t"))/len(end)'))
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

	# Separate Start from End
	outstart = []
	outend = []
	for row in rangeAT:
		pdElementstart = pd.DataFrame(item[0] for item in rangeAT[0])
		outstart.append(pdElementstart)
		pdElementend = pd.DataFrame(item[1] for item in rangeAT[1])
		outend.append(pdElementend)
	
	# Separate data frames by Start v End and Min v Max
	pdATCollectStartMin = pd.concat(outstart,axis=1).min(axis=1)
	pdATCollectEndMin = pd.concat(outend,axis=1).min(axis=1)
	pdATCollectStartMax = pd.concat(outstart,axis=1).max(axis=1)
	pdATCollectEndMax = pd.concat(outend,axis=1).max(axis=1)

	pdBins = pd.concat([pdATCollectStartMin,pdATCollectEndMin,pdATCollectStartMax,pdATCollectEndMax],axis=1)
	pdBins.columns=['StartMin','StopMin','StartMax','StopMax']
	return pdBins

def empericalATspread(element,size):
	totalSteps=[]
	for i in np.arange(1,size):
		pairStep = calculateAT(element,i)
		totalSteps.append(pairStep)
	return totalSteps

def main(rangeFeatures,fileName,binDir):
	directionBins,directionFeatures = evalN(rangeFeatures,fileName,binDir)
	return directionFeatures,directionBins

if __name__ == "__main__":
	main()