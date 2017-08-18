"""
Script to perform methylation analyses

Wren Saylor
July 5 2017

"""
import argparse
import pandas as pd
from ElementLibrary import get_bedtools_features
from ElementLibrary import convert_bedtool_to_panda
from ElementLibrary import save_bedtool_as_bedfile
from ElementLibrary import convert_panda_to_bed_format
from ElementLibrary import reverse_complement_dictionary
from ElementLibrary import get_just_fasta_sequence_for_feature
import GlobalVariables

# 	Old snippets from previous methylation processing method
# 	df.rename(columns=lambda x: x.replace('string','{0}_string'.format(filename)),inplace=True) # modify column names to be file specific
# 	df.set_index(keys='id',inplace=True,drop=True) # change index to column id
# 	df = df.iloc[::-1] # reverse order of index
# 	df = df.groupby(df.columns,axis=1).sum() # group by column name and sum those with same name
# 	dfstring = df.loc[:,df.columns.str.contains('string',case=False)] # search for column containing string and subset into new df
# 	df['group1'] = df.apply(lambda row:[i for i in row['group2'] if i in row['group3']],axis=1) # get a new column where two columns lists overlap
# 	df['string'] = df.apply(lambda row: [row['longstring'][i:i+2] for i in row['location']],axis=1) # get the local string for a location in the string
# 	df['dictionary'] = df.apply(lambda row: [dict(zip(row['group1'],row['group2']))],axis=1) # make a list a dictionary for another list in column
# 	List = df['list'].apply(pd.Series).stack().tolist() # make a list from a column of lists

# Threshold methylation data by coverage and percentage
def threshold_methylation_data(methFeatures):
	pdmethFeatures = convert_bedtool_to_panda(methFeatures)
	pdmethThresh = (pdmethFeatures[(pdmethFeatures.loc[:,3] >= GlobalVariables.methCovThresh) & (pdmethFeatures.loc[:,4] >= GlobalVariables.methPerThresh)])
	btmethThresh = convert_panda_to_bed_format(pdmethThresh)
	print 'Methylation coverage is being thresholded at {0} and percentage at {1}'.format(GlobalVariables.methCovThresh, GlobalVariables.methPerThresh)
	return btmethThresh

# Intersect regions from the methylation data with element regions
def intersect_methylation_data_by_element(rangeFeatures,methFeature):
	methSBoundary = methFeature.intersect(rangeFeatures[['chr','sBoundary','sEdge','id']].values.tolist(),wb=True,wa=True)
	if len(methSBoundary) != 0:
		pdmethSBoundary = convert_bedtool_to_panda(methSBoundary)
		pdmethSBoundary['int'] = 0
		pdmethSBoundary.columns = ['mchr','mstart','mstop','methCov','methPer','chr','sBoundary','sEdge','id','int']
		pdmethSBoundary['Context'] = pdmethSBoundary['mstop'] + 1
		pdmethSBoundary['BackContext'] = pdmethSBoundary['mstart'] -1
		pdmethSBoundary['Nuc'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethSBoundary[['mchr','BackContext','Context']].values.tolist()))
		pdmethSBoundary['NucContext'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethSBoundary[['mchr','mstop','Context']].values.tolist()))
		pdmethSBoundary['NucCytosine'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethSBoundary[['mchr','mstart','mstop']].values.tolist()))
		pdmethSBoundary['NucBackContext'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethSBoundary[['mchr','BackContext','mstart']].values.tolist()))
		pdmethSBoundary['methLoc'] = pdmethSBoundary['int'].astype(int)+(pdmethSBoundary['mstart'].astype(int)-pdmethSBoundary['sBoundary'].astype(int))
		outpdmethSBoundary = pdmethSBoundary[['chr','mstart','mstop','sBoundary','sEdge','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']]
		outpdmethSBoundary.columns = ['chr','methStart','methStop','eleStart','eleStop','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']
	else:
		outpdmethSBoundary = None

	methMiddle = methFeature.intersect(rangeFeatures[['chr','sCenter','eCenter','id']].values.tolist(),wb=True,wa=True)
	if len(methMiddle) != 0:
		pdmethFeature = convert_bedtool_to_panda(methMiddle)
		pdmethFeature['int'] = (((GlobalVariables.num - GlobalVariables.uce)/2) + GlobalVariables.inuce)
		pdmethFeature.columns = ['mchr','mstart','mstop','methCov','methPer','chr','sCenter','eCenter','id','int']
		pdmethFeature['Context'] = pdmethFeature['mstop'] + 1
		pdmethFeature['BackContext'] = pdmethFeature['mstart'] -1
		pdmethFeature['Nuc'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethFeature[['mchr','BackContext','Context']].values.tolist()))
		pdmethFeature['NucContext'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethFeature[['mchr','mstop','Context']].values.tolist()))
		pdmethFeature['NucCytosine'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethFeature[['mchr','mstart','mstop']].values.tolist()))
		pdmethFeature['NucBackContext'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethFeature[['mchr','BackContext','mstart']].values.tolist()))
		pdmethFeature['methLoc'] = pdmethFeature['int'].astype(int)+(pdmethFeature['mstart'].astype(int)-pdmethFeature['sCenter'].astype(int))
		outpdmethFeature = pdmethFeature[['chr','mstart','mstop','sCenter','eCenter','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']]
		outpdmethFeature.columns = ['chr','methStart','methStop','eleStart','eleStop','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']
	else:
		outpdmethFeature = None

	methEBoundary = methFeature.intersect(rangeFeatures[['chr','eEdge','eBoundary','id']].values.tolist(),wb=True,wa=True)
	if len(methEBoundary) != 0:
		pdmethEBoundary = convert_bedtool_to_panda(methEBoundary)
		pdmethEBoundary['int'] = GlobalVariables.num-1
		pdmethEBoundary.columns = ['mchr','mstart','mstop','methCov','methPer','chr','eEdge','eBoundary','id','int']
		pdmethEBoundary['Context'] = pdmethEBoundary['mstop'] + 1
		pdmethEBoundary['BackContext'] = pdmethEBoundary['mstart'] -1
		pdmethEBoundary['Nuc'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethEBoundary[['mchr','BackContext','Context']].values.tolist()))
		pdmethEBoundary['NucContext'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethEBoundary[['mchr','mstop','Context']].values.tolist()))
		pdmethEBoundary['NucCytosine'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethEBoundary[['mchr','mstart','mstop']].values.tolist()))
		pdmethEBoundary['NucBackContext'] = get_just_fasta_sequence_for_feature(get_bedtools_features(pdmethEBoundary[['mchr','BackContext','mstart']].values.tolist()))
		pdmethEBoundary['methLoc'] = pdmethEBoundary['int'].astype(int)-(pdmethEBoundary['eBoundary'].astype(int)-pdmethEBoundary['mstop'].astype(int))
		outpdmethEBoundary = pdmethEBoundary[['id','methPer','methLoc','methCov']]
		outpdmethEBoundary = pdmethEBoundary[['chr','mstart','mstop','eEdge','eBoundary','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']]
		outpdmethEBoundary.columns = ['chr','methStart','methStop','eleStart','eleStop','int','id','methPer','methLoc','methCov','Nuc','NucContext','NucCytosine','NucBackContext']
	else:
		outpdmethEBoundary = None

	methList = [outpdmethSBoundary,outpdmethFeature,outpdmethEBoundary]
	concatMeth = pd.concat(methList)
	sortMeth = concatMeth.sort_values(['methLoc'],ascending=True)
	return sortMeth

# Run the analysis to extract percentage, frequency, coverage, location, context, and direction
def collect_methylation_data_by_element(rangeFeatures):
	outMeth = []
	for methName in GlobalVariables.mFiles:
		methFeatures = get_bedtools_features(methName)
		pdmethThresh = threshold_methylation_data(methFeatures)
		methPosition = intersect_methylation_data_by_element(rangeFeatures,pdmethThresh)
		methPosition['tissue'] = methName.replace('.bed','')
		stringDF = rangeFeatures[['id','combineString']]
		methMerge = pd.merge(methPosition,stringDF,how='left',on='id')
		methMerge['methLocBEnd'] = methMerge['methLoc'] - 1
		methMerge['methLocCEnd'] = methMerge['methLoc'] + 1
		methMerge['methLocEnd'] = methMerge['methLoc'] + 2
		methMerge['Cytosine'] = methMerge.apply(lambda row: row['combineString'][row['methLoc']:row['methLocCEnd']],axis=1)
		methMerge['Context'] = methMerge.apply(lambda row: row['combineString'][row['methLocCEnd']:row['methLocEnd']],axis=1)
		methMerge['BackContext'] = methMerge.apply(lambda row: row['combineString'][row['methLocBEnd']:row['methLoc']],axis=1)
		methMerge['ContextCheck'] = methMerge.apply(lambda row: row['combineString'][row['methLocBEnd']:row['methLocEnd']],axis=1)
		methMerge['methFreq'] = methMerge.groupby(['methLoc','Cytosine'])['methLoc'].transform('count')
		methMerge['Nuc'] = methMerge['Nuc'].str.upper()
		methMerge['sameSeq'] = methMerge['Nuc'] == methMerge['ContextCheck']
		
		# If the nucleotide in the cytosine column is 'G', make the context the other direction (reverse complement later, in graphing, in order to differentiate between strands)
		methMerge.loc[methMerge['Cytosine'] == 'G', 'Context'] = methMerge['BackContext']
		
		# sameSeq might be 'False' if 1) the c is at the end border for the downstream boundary, 2) the sequence bridges the sequence split for the upstream boundary
		falseMeth = (methMerge[methMerge['sameSeq'] == False])
		
		# Conditionally update contexts where False for matches between sequence and methylation nucleotides- c get context, and g gets backcontext
		methMerge.loc[methMerge['sameSeq'] == False,'Cytosine'] = methMerge['NucCytosine']
		methMerge.loc[(methMerge['sameSeq'] == False) & (methMerge['NucCytosine'] == 'C'),'Context'] = methMerge['NucContext']
		methMerge.loc[(methMerge['sameSeq'] == False) & (methMerge['NucCytosine'] == 'G'),'Context'] = methMerge['NucBackContext']
		
		print 'There are {0} instances at {1} where methylation context did not match between methylation bedfile and sequence in {2}'.format(len(falseMeth.index),falseMeth['methLoc'].tolist(),methName)
		subMeth = methMerge[['id','methLoc','methPer','methCov','methFreq','Cytosine','Context','tissue']]
		outMeth.append(subMeth)
	print 'Discrepancy of context are acceptable at the end of the sequence, or at the split between the middle and boundary, the context extracted from the original methylation file will be used instead'
	pdMeth = pd.concat(outMeth)
	return pdMeth

def main(rangeFeatures):
	print 'Running MethlationLibrary'
	pdMeth = collect_methylation_data_by_element(rangeFeatures)
	return pdMeth

if __name__ == "__main__":
	main()