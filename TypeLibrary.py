"""
Script to separate by type

Wren Saylor
July 5 2017
"""
import argparse
import pandas as pd
from MethylationLibrary import compactMeth
from FangsLibrary import compactWindow
import GlobalVariables

# do all the analysis for each type
def perType(boolType,fileName):
	typeWindow, typeNames = compactWindow(boolType['combineString'],boolType['id'])
	if any(x in GlobalVariables.graphs for x in ['methylation','cluster']):
		pdMeth = compactMeth(boolType)
	else: 
		pdMeth = None
	return pdMeth,typeWindow,typeNames

def main(boolType,fileName):
	print 'Running TypeLibrary'
	pdMeth,typeWindow,typeNames = perType(boolType,fileName)
	print 'Completed group sorting for {0} items'.format(len(boolType.index))
	return pdMeth,typeWindow,typeNames

if __name__ == "__main__":
	main() 