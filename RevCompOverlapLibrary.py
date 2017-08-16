"""
Script to do the reverse complement sorting for the aligned overlap regions

Wren Saylor
August 10 2017
"""

import argparse
import pandas as pd
from RevCompLibrary import methDirection
from RevCompLibrary import slideDirection
import GlobalVariables

# RCsort and return methylation and sliding window computations
def dirOverlapsLine(posStr,negStr):
	compWindow, compNames = slideDirection(negStr,posStr)
	if any(x in GlobalVariables.graphs for x in ['methylation','cluster']):
		groupMeth = methDirection(negStr,posStr)
	else: 
		groupMeth = None
	return groupMeth,compWindow,compNames

def main(posStr,negStr):
	print 'Running RevCompOverlapLibrary'
	groupMeth,compWindow,compNames = dirOverlapsLine(posStr,negStr)
	return groupMeth,compWindow,compNames

if __name__ == "__main__":
	main()