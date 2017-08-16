"""
Script to perform sliding window to return content for nucleotides

Wren Saylor
"""

import argparse
import pandas as pd
import GlobalVariables

# run the sliding window for each nucleotide string
def compactWindow(features,label):
	outCollect = []
	for element,id in zip(features,label):
		outElement = {id: []}
		outList = {key:[] for key in GlobalVariables.nucList}
		n = GlobalVariables.num
		s = 1 # size to jump for sliding window
		start, end = 0, GlobalVariables.window
		while end < n:
			current = element[start:end]
			for key in GlobalVariables.nucList:
				percentage = float(100*current.count(key)/len(current))
				outList[key].append(percentage)
			start, end = start + s, end + s
		outElement[id].append(outList)
		outCollect.append(outElement)
	outFlatten = flattenWindow(outCollect)
	outDataFrame, names = dfWindow(outFlatten)
	print 'Retrieved sliding window data for nucleotides strings {0}'.format(names)
	return outDataFrame,names

# convert to a data frame with a list for each element x nucleotide string
def flattenWindow(outCollect):
	outFlatten = pd.DataFrame()
	outIndex = []
	for d in outCollect:
		for k,v in d.items():
			outFlatten = outFlatten.append(v,ignore_index=True)
			outIndex.append(k)
	outFlatten.index = outIndex
	return outFlatten

# turn each list of element x nucleotide string into a separate df, within a larger df
def dfWindow(outDataFrame):
	names = outDataFrame.columns.tolist()
	collectNucDF = []
	for nuc in names:
		nuc = outDataFrame[[nuc]]
		nuc.columns = ['temp']
		split = pd.DataFrame(nuc['temp'].values.tolist(),index=nuc.index)
		collectNucDF.append(split)
	return collectNucDF,names

def main(features,label):
	print 'Running FangsLibrary'
	outDataFrame, names = compactWindow(features,label)
	return outDataFrame, names

if __name__ == "__main__":
	main()