"""
Script to do the reverse complement sorting for the aligned overlap regions

Wren Saylor
August 10 2017
"""

import argparse
import pandas as pd
from RevCompLibrary import methDirection
from RevCompLibrary import slideDirection

# RCsort and return methylation and sliding window computations
def dirOverlapsLine(posStr,negStr,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs):
	compWindow, compNames = slideDirection(negStr,posStr,num,uce,inuce,window,nucLine)
	if any(x in graphs for x in ['methylation','cluster']):
		groupMeth = methDirection(negStr,posStr,mFiles,num,uce,inuce,methCovThresh,methPerThresh,faGenome)
	else: 
		groupMeth = None
	return groupMeth,compWindow,compNames

def main(posStr,negStr,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs):
	print 'Running RevCompOverlapLibrary'
	groupMeth,compWindow,compNames = dirOverlapsLine(posStr,negStr,mFiles,num,uce,inuce,window,methCovThresh,methPerThresh,nucLine,faGenome,graphs)
# 	print 'Completed reverse complements sorting for {0}  items and {1} neg items'.format(len(posStr,negStr.index))
	return groupMeth,compWindow,compNames

if __name__ == "__main__":
	main()