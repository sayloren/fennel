"""
Script to separate by directionality

Wren Saylor
July 5 2017
"""
import argparse
import numpy as np
import pandas as pd
import GlobalVariables

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
	rangeFeatures['compareBoundaries'] = rangeFeatures.apply(lambda row: (compare_boundaries_size_n(row['feature'],GlobalVariables.binDir)),axis=1)
	rangeFeatures['compareBoundariesRange'] = rangeFeatures.apply(lambda row: (emperical_boundaries_for_increasing_bins_symbol(row['feature'],GlobalVariables.binDir)),axis=1)
	rangeFeatures['equalBoundariesCount'] = rangeFeatures.apply(lambda row: row['compareBoundariesRange'].count('='),axis=1)
	rangeFeatures['plusBoundariesCount'] = rangeFeatures.apply(lambda row: row['compareBoundariesRange'].count('+'),axis=1)
	rangeFeatures['minusBoundariesCount'] = rangeFeatures.apply(lambda row: row['compareBoundariesRange'].count('-'),axis=1)
# 	print rangeFeatures[['plusBoundariesCount','minusBoundariesCount','equalBoundariesCount']]
	
	compareEnds = pd.DataFrame(rangeFeatures[['chr','start','end','compareBoundaries']])
	print 'Sorting the element boundaries by bin size {0}'.format(GlobalVariables.binDir)
	rangeBins = collect_emperical_boundary_comparisons(rangeFeatures)
	return rangeFeatures,rangeBins

# Get the actual spread of = in the data over increasing bin size
def collect_emperical_boundary_comparisons(rangeFeatures):
	# Perform min,max collection
	rangeAT = rangeFeatures.apply(lambda row: (emperical_boundaries_for_increasing_bins_percentage(row['feature'],GlobalVariables.binDir)),axis=1)
	pdrangeAT = pd.DataFrame(rangeAT.values.tolist())
	
	# Get direction collection
	equalAT = rangeFeatures.apply(lambda row: (emperical_boundaries_for_increasing_bins_symbol(row['feature'],GlobalVariables.binDir)),axis=1)
	pdequalAT = pd.DataFrame(equalAT.values.tolist())
	countAT = pdequalAT.apply(pd.value_counts)
	equalCounts = countAT.loc['=']/len(pdequalAT.index)
	
	# Separate Start from End for min and max
	splitAT = pd.concat(dict([(row[0],row[1].apply(lambda y: pd.Series(y))) for row in pdrangeAT.iterrows()]),axis=1)
	
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

	pdBins = pd.concat([Min,Max,equalCounts],axis=1)
	pdBins.columns=['Min','Max','Equal']
	return pdBins

# Run calculate_nucleotides_at for increasing bin sizes
def emperical_boundaries_for_increasing_bins_percentage(element,size):
	totalSteps=[]
	for i in np.arange(1,size*2):
		pairStep = calculate_nucleotides_at(element,i)
		totalSteps.append(pairStep)
	return totalSteps

# Run compare_boundaries_size_n over increasing bin sizes
def emperical_boundaries_for_increasing_bins_symbol(element,size):
	totalSteps=[]
	for i in np.arange(1,size*2):
		perSize = calculate_nucleotides_at(element,i)
		# give + - = depending on which side has larger AT content
		if perSize[0] > perSize[1]: outList = '+'
		if perSize[1] > perSize[0]: outList = '-'
		if perSize[1] == perSize[0]: outList = '='
		totalSteps.append(outList)
	return totalSteps

def main(rangeFeatures,fileName):
	print 'Running DirectionLibrary'
	directionFeatures,directionBins = evaluate_boundaries_size_n(rangeFeatures,fileName)
	return directionFeatures,directionBins

if __name__ == "__main__":
	main()