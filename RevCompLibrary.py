"""
Script to perform RC sorting

Wren Saylor
July 5 2017

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
import pandas as pd
from FangsLibrary import run_sliding_window_for_each_nucleotide_string
from ElementLibrary import get_bedtools_features
from MethylationLibrary import collect_methylation_data_by_element
import GlobalVariables

# Methylation RCsorting
def sort_methylation_by_directionality(negStr,posStr):

	posMeth = collect_methylation_data_by_element(posStr)
	negMeth = collect_methylation_data_by_element(negStr)
	
	# Zip reversed range to make a dictionary for replacing the location of the neg methylation
	originalRange = range(0,GlobalVariables.num)
	reverseRange = originalRange[::-1]
	rangeDict = dict(zip(originalRange,reverseRange))
	
	# Zip reverse complement sequence for replacing the nucleotides for neg methylation
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}

	# Convert neg Meth df
	negMeth['methLocNew'] = negMeth.methLoc.map(rangeDict)
	negMeth['CytosineNew'] = negMeth.Cytosine.map(seqDict)
	negMeth['ContextNew'] = negMeth.Context.map(seqDict)
	negMethNew = negMeth[['id','methLocNew','methPer','methCov','methFreq','CytosineNew','ContextNew','tissue']]
	negMethNew.columns = ['id','methLoc','methPer','methCov','methFreq','Cytosine','Context','tissue']
	
	# Concat pos and revised neg meth dfs
	frames = [posMeth,negMethNew]
	catMerge = pd.concat(frames)
	
	# Update Frequencey count column
	catMerge['methFreqNew'] = catMerge.groupby(['methLoc','tissue','Cytosine'])['methLoc'].transform('count')
	outMerge = catMerge[['id','methLoc','methPer','methCov','methFreqNew','Cytosine','Context','tissue']]
	outMerge.columns = ['id','methLoc','methPer','methCov','methFreq','Cytosine','Context','tissue']

	return outMerge

# Sliding window RCsorting
def sort_sliding_window_by_directionality(negStr,posStr):
	negDF, negNames = run_sliding_window_for_each_nucleotide_string(negStr['reverseComplement'],negStr['id'])
	posDF, posNames = run_sliding_window_for_each_nucleotide_string(posStr['combineString'],posStr['id'])
	compWindow = []
	for x, y in zip(negDF, posDF):
		tempCat = pd.concat([x,y],axis=1)
		tempGroup = tempCat.groupby(tempCat.columns,axis=1).sum()
		compWindow.append(tempGroup)
	return compWindow, negNames

# Separate on plus and minus orientation, RCsort and return methylation and sliding window computations
def sort_elements_by_directionality(directionFeatures):
	negStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '-')])
	posStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '+')])
	compWindow, compNames = sort_sliding_window_by_directionality(negStr,posStr)
	if any(x in GlobalVariables.graphs for x in ['methylation','cluster','methextend']):
		groupMeth = sort_methylation_by_directionality(negStr,posStr)
	else: 
		groupMeth = None
	return groupMeth,compWindow,compNames

def main(directionFeatures):
	groupMeth,compWindow,compNames = sort_elements_by_directionality(directionFeatures)
	print 'Completed reverse complement sorting for {0} items, with {1} bin sorting'.format(len(directionFeatures.index),GlobalVariables.binDir)
	return groupMeth,compWindow,compNames

if __name__ == "__main__":
	main()
