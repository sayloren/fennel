"""
Script to perform sliding window to return content for nucleotides

Wren Saylor

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
import GlobalVariables

# run the sliding window for each nucleotide string
def run_sliding_window_for_each_nucleotide_string(features,label):
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
	outFlatten = flatten_data_from_sliding_window(outCollect)
	outDataFrame, names = convert_sliding_window_to_dataframe(outFlatten)
	print 'Retrieved sliding window data for nucleotides strings {0}'.format(names)
	return outDataFrame,names

# convert to a data frame with a list for each element x nucleotide string
def flatten_data_from_sliding_window(outCollect):
	outFlatten = pd.DataFrame()
	outIndex = []
	for d in outCollect:
		for k,v in d.items():
			outFlatten = outFlatten.append(v,ignore_index=True)
			outIndex.append(k)
	outFlatten.index = outIndex
	return outFlatten

# turn each list of element x nucleotide string into a separate df, within a larger df
def convert_sliding_window_to_dataframe(outDataFrame):
	names = outDataFrame.columns.tolist()
	collectNucDF = []
	for nuc in names:
		nuc = outDataFrame[[nuc]]
		nuc.columns = ['temp']
		split = pd.DataFrame(nuc['temp'].values.tolist(),index=nuc.index)
		collectNucDF.append(split)
	return collectNucDF,names

def main(features,label):
	outDataFrame, names = run_sliding_window_for_each_nucleotide_string(features,label)
	return outDataFrame, names

if __name__ == "__main__":
	main()
