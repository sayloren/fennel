"""
Script to do the reverse complement sorting for the aligned overlap regions

Wren Saylor
August 10 2017
"""

import argparse
import pandas as pd
from RevCompLibrary import sort_methylation_by_directionality
from RevCompLibrary import sort_sliding_window_by_directionality
import GlobalVariables

# RCsort and return methylation and sliding window computations
def run_directionality_for_overlaps(posStr,negStr):
	compWindow, compNames = sort_sliding_window_by_directionality(negStr,posStr)
	if any(x in GlobalVariables.graphs for x in ['methylation','cluster']):
		groupMeth = sort_methylation_by_directionality(negStr,posStr)
	else: 
		groupMeth = None
	return groupMeth,compWindow,compNames

def main(posStr,negStr):
	print 'Running RevCompOverlapLibrary'
	groupMeth,compWindow,compNames = run_directionality_for_overlaps(posStr,negStr)
	return groupMeth,compWindow,compNames

if __name__ == "__main__":
	main()