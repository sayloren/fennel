"""
Script to align exons

Wren Saylor
August 4 2017

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
from collections import defaultdict
import pandas as pd
import numpy as np
import pybedtools as pbt
from ElementLibrary import reverse_complement_dictionary
from ElementLibrary import get_bedtools_features
from ElementLibrary import convert_bedtool_to_panda
from ElementLibrary import get_just_fasta_sequence_for_feature
import GlobalVariables

def save_panda_data_frame(pdData, strFilename):
	pdData.to_csv(strFilename, sep='\t', header=False, index=False)

def make_data_frame_for_id_type_label(df,fileName):
	df['crossboundary'] = np.where((df['startdifference'] >= 0) & (df['enddifference']  >= 0),'No','Yes') # if elemenet complete within exon, label no, else yes
	outcross = df[['id','crossboundary']]
	save_panda_data_frame(outcross, 'ElementCrossBoundary_{0}.txt'.format(fileName))

# Intersect the exons with the uces
def intersect_overlaping_regions_by_element(df,overlap,fileName):
	OverlapBoundary = overlap.intersect(df[['chr','start','end','id']].values.tolist(),wo=True)
	pdOverlapBoundary = convert_bedtool_to_panda(OverlapBoundary)
	pdOverlapBoundary.columns = ['ochr','ostart','oend','chr','start','end','id','overlap']#,'exondir'
	pdOverlapBoundary['elementsize'] = pdOverlapBoundary['end'] - pdOverlapBoundary['start']
	pdOverlapBoundary['startdifference'] = pdOverlapBoundary['start'] - pdOverlapBoundary['ostart']# if this is a negative number, the element starts before the exon
	pdOverlapBoundary['enddifference'] = pdOverlapBoundary['oend'] - pdOverlapBoundary['end']# if this is a negative number, the element ends after the exon
	overlapTable = make_table_for_overlap_classification(pdOverlapBoundary,fileName)
	return pdOverlapBoundary, overlapTable

# Make a table for how many of each type of overlapping region
def make_table_for_overlap_classification(df,fileName):
	make_data_frame_for_id_type_label(df,fileName)
	complete = df[(df['startdifference'] >= 0) & (df['enddifference']  >= 0)].count() # how many elements are completely inside the overlaps
	start = df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)].count() # how many elements overlap at the upstream boundary
	end = df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)].count() # how many elements overlap at the downstream boundary
	interior = df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)].count() # how many elements contain an overlaps
	crossboundary = start + end
	multi = df[df.duplicated(subset='id',keep='first')].count() # number of elements that fit into two categories
	multiid = df[df.duplicated(subset='id',keep=False)].groupby('id').min()
	print 'Elements {0} where represented in two groups'.format(multiid.index)
	total = df.count() - multi # get the total number of elements intersecting overlaps - minus those that overlap two overlaps
	frames = [complete,crossboundary,interior,total]
	overlapTable = pd.concat(frames,axis=1)
	overlapTable.columns = ['embedded','crossboundary','contains','total']
	overlapTable = overlapTable.head(1)
	overlapTable.index = ['overlaps']
	return overlapTable

# Get the range for each separate classification of overlapping regions
def classify_overlap_regions_by_type_overlap(df):
	centerElement = GlobalVariables.num/2
	flankSize = (GlobalVariables.num - GlobalVariables.uce)/2
	inregion = GlobalVariables.uce-(GlobalVariables.inuce*2)
	# Overlap regions that are completely within an overlaps
	if len(df[(df['startdifference'] >= 0) & (df['enddifference'] >= 0)]) != 0: # elements are completely inside the overlaps
		features = df[(df['startdifference'] >= 0) & (df['enddifference'] >= 0)]
		features.reset_index(drop=True,inplace=True)
		features['middle'] = features[['start','end']].mean(axis=1).astype(int)
		features['sCenter'] = features.loc[:,'middle'] - (inregion/2)
		features['eCenter'] = features.loc[:,'middle'] + (inregion/2)
		features['sEdge'] = features.loc[:,'start'] + GlobalVariables.inuce
		features['eEdge'] = features.loc[:,'end'] - GlobalVariables.inuce
		features['sBoundary'] = features.loc[:,'start'] - flankSize
		features['eBoundary'] = features.loc[:,'end'] + flankSize
		features['sBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','sBoundary','sEdge']].values.tolist()))
		features['MiddleSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','sCenter','eCenter']].values.tolist()))
		features['eBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','eEdge','eBoundary']].values.tolist()))
		features['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','start','end']].values.tolist()))
		features['combineString'] = features['sBoundarySeq'].astype(str) +features['MiddleSeq'].astype(str) + features['eBoundarySeq'].astype(str)
		features['combineString'] = features['combineString'].str.upper()
		features['size'] = features['combineString'].str.len() # check length of string
		features['reverseComplement'] = features.apply(lambda row: reverse_complement_dictionary(row['combineString']),axis=1)
		completeelement = features[['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','feature','id']]
		completeelement.columns = ['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','feature','id']
	else:
		completeelement = None
	# Exonic regions that contain an overlaps - was going to align by interior overlaps, buuuut, its tricky
	if len(df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)]) != 0: # elements contain an overlaps
		features = df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)]
		interiorsize = df['overlap'].min() # smallest overlaps, to be able to compare the exon boundaries crossover
		df['difference'] = df['elementsize'] - df['overlap']
		elementsize = df['difference'].min()
		features = df[(df['startdifference'] <= -1) & (df['enddifference'] <= -1)]
		features.reset_index(drop=True,inplace=True)
		features['middle'] = features[['start','end']].mean(axis=1).astype(int)
		features['sCenter'] = features['middle'] - (inregion/2)
		features['eCenter'] = features['middle'] + (inregion/2)
		features['sEdge'] = features['start'] + GlobalVariables.inuce
		features['eEdge'] = features['end'] - GlobalVariables.inuce
		features['sBoundary'] = features['start'] - flankSize
		features['eBoundary'] = features['end'] + flankSize
		features['sBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','sBoundary','start']].values.tolist()))
		features['sEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','start','sEdge']].values.tolist()))
		features['MiddleSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','sCenter','eCenter']].values.tolist()))
		features['eEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','eEdge','end',]].values.tolist()))
		features['eBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','end','eBoundary']].values.tolist()))
		features['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','start','end']].values.tolist()))
		features['combineString'] = features['sBoundarySeq'].astype(str) + features['sEdgeSeq'].astype(str) + features['MiddleSeq'].astype(str) + features['eEdgeSeq'].astype(str) + features['eBoundarySeq'].astype(str)
		features['combineString'] = features['combineString'].str.upper()
		features['size'] = features['combineString'].str.len() # check length of string
		features['reverseComplement'] = features.apply(lambda row: reverse_complement_dictionary(row['combineString']),axis=1)
		interiorelement = features[['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','feature','id']]
		interiorelement.columns = ['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','feature','id']
	else:
		interiorelement = None
	# Exonic cross boundary shift (still need to separate by intergenic/intronic)
	if len(df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)]) != 0: # elements overlap at the upstream boundary
		features = df[(df['startdifference'] <= -1) & (df['enddifference'] >= 0)]
		features.reset_index(drop=True,inplace=True)
		features['newstart'] = features['oend']
		features['newend'] = features['oend'] + GlobalVariables.uce
		features['sEdge'] = features['newstart'] + GlobalVariables.inuce
		features['eEdge'] = features['newend'] - GlobalVariables.inuce
		features['sCenter'] = features['newstart'] + (inregion/2)
		features['eCenter'] = features['newend'] - (inregion/2)
		features['sBoundary'] = features['newstart'] - flankSize
		features['eBoundary'] = features['newend'] + flankSize
		features['combineString'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','sBoundary','eBoundary']].values.tolist()))
		features['combineString'] = features['combineString'].str.upper()
		features['size'] = features['combineString'].str.len() # check length of string
		features['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','newstart','newend']].values.tolist()))
		features['reverse_complement_dictionary'] = features.apply(lambda row: reverse_complement_dictionary(row['combineString']),axis=1)
		startcrossboundary = features[['chr','sBoundary','newstart','sEdge','sCenter','eCenter','eEdge','newend','eBoundary','combineString','reverse_complement_dictionary','feature','id']]
		startcrossboundary.columns = ['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverse_complement_dictionary','feature','id']
	else:
		startcrossboundary = None
	if len(df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)]) != 0: # elements overlap at the downstream boundary
		features = df[(df['startdifference'] >= 0) & (df['enddifference'] <= -1)]
		features.reset_index(drop=True,inplace=True)
		features['newstart'] = features['ostart'] - GlobalVariables.uce
		features['newend'] = features['ostart']
		features['sEdge'] = features['newstart'] + GlobalVariables.inuce
		features['eEdge'] = features['newend'] - GlobalVariables.inuce
		features['sCenter'] = features['newstart'] + (inregion/2)
		features['eCenter'] = features['newend'] - (inregion/2)
		features['sBoundary'] = features['newstart'] - flankSize
		features['eBoundary'] = features['newend'] + flankSize
		features['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','newstart','newend']].values.tolist()))
		features['combineString'] = get_just_fasta_sequence_for_feature(get_bedtools_features(features[['chr','sBoundary','eBoundary']].values.tolist()))
		features['combineString'] = features['combineString'].str.upper()
		features['reverseComplement'] = features.apply(lambda row: reverse_complement_dictionary(row['combineString']),axis=1)
		features['size'] = features['combineString'].str.len() # check length of string
		endcrossboundary = features[['chr','sBoundary','newstart','sEdge','sCenter','eCenter','eEdge','newend','eBoundary','combineString','reverseComplement','feature','id']]
		endcrossboundary.columns = ['chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary','combineString','reverseComplement','feature','id']
	else:
		endcrossboundary = None
	return startcrossboundary,endcrossboundary,completeelement,interiorelement

def main(rangeFeatures,fileName):
	print 'Running OverlapLibrary'
	overFeature = get_bedtools_features(GlobalVariables.Overlapregions)
	overlapFeatures, overlapTable = intersect_overlaping_regions_by_element(rangeFeatures,overFeature,fileName)
	startcrossboundary,endcrossboundary,completeelement,interiorelement = classify_overlap_regions_by_type_overlap(overlapFeatures)
	return startcrossboundary,endcrossboundary,completeelement,interiorelement,overlapTable

if __name__ == "__main__":
	main()
