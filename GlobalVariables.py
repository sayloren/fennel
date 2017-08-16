"""
Script to set global variables

Wren Saylor
August 15 2017
"""
import argparse

def setGlobals(args):
	# Integer parameters
	global num
	global uce
	global inuce
	global window
	global binDir
	global methCovThresh
	global methPerThresh
	global methFlank
	global halfwindow
	global fillX
	num = args.total
	uce = args.element
	inuce = args.inset
	window = args.window
	binDir = args.bin
	methCovThresh = args.thresholdcoverage
	methPerThresh = args.thresholdpercentage
	methFlank = args.methylationflank
	halfwindow = ((window/2)+1)
	fillX = range(0,(num-window))

	# Element, random regions and methylation files
	global eFiles
	global mFiles
	global rFiles
	eFiles = [line.strip() for line in args.efile]
	mFiles = [line.strip() for line in args.mfile]
	rFiles = [line.strip() for line in args.rfile]

	# Genome files from UCSC
	global sizeGenome
	global faGenome
	global Overlapregions
	sizeGenome = args.genome
	faGenome = args.fasta
	Overlapregions = args.overlapingelements

	# Lists with the types and directions to use
	global typeList
	global dirList
	global nucList
	typeList = args.elementype
	dirList = args.elementdirection
	nucList = args.nucleotideline

	# Reverse complement argument
	global revCom
	global overlapInset
	revCom = args.reversecomplement
	overlapInset = args.elementalign

	# Which plots to run
	global graphs
	graphs = args.plots

	# A string to add to the out file name in case you want to set up runs and let be
	global stringName
	stringName = args.stringname
	
	# Locations for plotting with sliding window
	global plotLineLocationOne # upstream element boundary
	global plotLineLocationTwo # downstream element boundary
	global plotLineLocationThree # upstream element inset
	global plotLineLocationFour # downstream element inset
	plotLineLocationOne = (((num-uce)/2)+(inuce-halfwindow))
	plotLineLocationTwo = (((num-uce)/2)+(uce-inuce-halfwindow))
	plotLineLocationThree = (((num-uce)/2)-halfwindow)
	plotLineLocationFour = (((num-uce)/2)+uce-halfwindow)
	
	# Locations for plotting without the sliding window
	global plotLineLocationOneFull # upstream element boundary
	global plotLineLocationTwoFull # downstream element boundary
	global plotLineLocationThreeFull # upstream element inset
	global plotLineLocationFourFull # downstream element inset
	plotLineLocationOneFull = (((num-uce)/2)+inuce)
	plotLineLocationTwoFull = (((num-uce)/2)+(uce-inuce))
	plotLineLocationThreeFull = ((num-uce)/2)
	plotLineLocationFourFull = (((num-uce)/2)+uce)

def main(args):
	setGlobals(args)

if __name__ == "__main__":
	main()