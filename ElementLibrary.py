"""
Script to process the elements you what to analyze

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
from collections import defaultdict
import pandas as pd
import pybedtools as pbt
import GlobalVariables

# get bt features
def get_bedtools_features(strFileName):
	btFeatures = pbt.BedTool(strFileName)
	return btFeatures

# get the correct range for fang evaluation
def collect_coordinates_for_element_positions(btFeatures,fileName):#btFeatures,fileName,num,uce,inuce
	flankSize = (GlobalVariables.num - GlobalVariables.uce)/2
	inregion = GlobalVariables.uce-(GlobalVariables.inuce*2)
	midFeatures = pd.read_table(btFeatures.fn, header=None)
	midFeatures['middle'] = midFeatures.loc[:,1:2].mean(axis=1).astype(int)
	midFeatures['start'] = midFeatures.loc[:,1]
	midFeatures['end'] = midFeatures.loc[:,2] 
	midFeatures['sCenter'] = midFeatures['middle'] - (inregion/2)
	midFeatures['eCenter'] = midFeatures['middle'] + (inregion/2)
	midFeatures['sEdge'] = midFeatures.loc[:,1] + GlobalVariables.inuce
	midFeatures['eEdge'] = midFeatures.loc[:,2] - GlobalVariables.inuce
	midFeatures['sBoundary'] = midFeatures.loc[:,1] - flankSize
	midFeatures['eBoundary'] = midFeatures.loc[:,2] + flankSize
	midFeatures['type'] = midFeatures.loc[:,4]
	midFeatures['id'] = midFeatures.loc[:,3]
	midFeatures['chr'] = midFeatures.loc[:,0]
	midFeatures['size'] = midFeatures.loc[:,2].astype(int)-midFeatures.loc[:,1].astype(int)
	rangeFeatures = midFeatures[['type','id','size','chr','sBoundary','start','sEdge','sCenter','eCenter','eEdge','end','eBoundary']]
	print 'Getting the coordinates for the area to examine, {0} from the middle of the element, {1} inset from the edges'.format((midFeatures['eCenter']-midFeatures['sCenter']).mean(),(midFeatures['sEdge']-midFeatures['start']).mean())
	return rangeFeatures

# get the strings for sliding window regions
def get_fasta_for_element_coordinates(rangeFeatures):#rangeFeatures,faGenome
	rangeFeatures['sBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sBoundary','start']].values.tolist()))
	rangeFeatures['sEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','sEdge']].values.tolist()))
	rangeFeatures['MiddleSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','sCenter','eCenter']].values.tolist()))
	rangeFeatures['eEdgeSeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','eEdge','end',]].values.tolist()))
	rangeFeatures['eBoundarySeq'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','end','eBoundary']].values.tolist()))
	rangeFeatures['feature'] = get_just_fasta_sequence_for_feature(get_bedtools_features(rangeFeatures[['chr','start','end']].values.tolist()))
	rangeFeatures['combineString'] = rangeFeatures['sBoundarySeq'].astype(str) + rangeFeatures['sEdgeSeq'].astype(str) + rangeFeatures['MiddleSeq'].astype(str) + rangeFeatures['eEdgeSeq'].astype(str) + rangeFeatures['eBoundarySeq'].astype(str)
	rangeFeatures['combineString'] = rangeFeatures['combineString'].str.upper()
	rangeFeatures['feature'] = rangeFeatures['feature'].str.upper()
	rangeFeatures['reverseComplement'] = rangeFeatures.apply(lambda row: reverse_complement_dictionary(row['combineString']),axis=1)
	return rangeFeatures

# Get the percentage AT in the element
def percentage_at_for_element(region,fileName):
	collectAT = []
	for r in region:
		collectAT.append(eval('100*float(r.count("A") + r.count("a") + r.count("T") + r.count("t"))/len(r)'))
	pdAT = pd.DataFrame(collectAT)
	print 'Mean AT content for all {0} elements in {1} is {2} Percent'.format(len(region.index),fileName,pdAT.mean())

# get the reverse complement
def reverse_complement_dictionary(sequence):
	seqDict = {'A':'T','T':'A','C':'G','G':'C','N':'N'}
	return "".join([seqDict[base] for base in reversed(sequence)])

# save file from bedtool
def save_bedtool_as_bedfile(btObject,strFilename):
	btObject.saveas(strFilename)

# convert bedtool to panda
# If there is nothing in the btobject, will it read the data from the previous itteration!?
def convert_bedtool_to_panda(btobject):
	save_bedtool_as_bedfile(btobject,'temp.bed')
	pdObject = pd.read_table(btobject.fn, header=None)
	return pdObject

# convert panda to bedtool
def convert_panda_to_bed_format(panda):
	arArFeatures = panda.values.tolist()
	btoutFeatures = get_bedtools_features(arArFeatures)
	return btoutFeatures

# get fasta strings for each desired region
def get_fast_for_type(btFeatures,fileName):
	saveAll = 'Seq_result_for_all_{0}.txt'.format(fileName)
	seqAll = btFeatures.sequence(fi=GlobalVariables.faGenome,fo=saveAll)
	saveExonic = 'Seq_results_for_exonic_{0}.txt'.format(fileName)
	seqExonic = btFeatures.sequence(fi=GlobalVariables.faGenome,fo=saveExonic).filter(lambda x: x[name] == 'exonic')
	saveIntronic = 'Seq_results_for_intronic_{0}.txt'.format(fileName)
	seqIntronic = btFeatures.sequence(fi=GlobalVariables.faGenome,fo=saveIntronic).filter(lambda x: x[name] == 'intronic')
	saveIntergenic = 'Seq_results_for_intergenic_{0}.txt'.format(fileName)
	seqIntergenic = btFeatures.sequence(fi=vfaGenome,fo=saveIntergenic).filter(lambda x: x[name] == 'intergenic')
	return saveAll, saveExonic, saveIntronic, saveIntergenic

# used in get_fasta_for_element_coordinates to extract just the fasta strings
def get_just_fasta_sequence_for_feature(inFeature):
	seqFeature = inFeature.sequence(fi=GlobalVariables.faGenome)
	outFeature = pd.read_table(seqFeature.seqfn)
	outSequence = outFeature[::2]
	outSequence = outSequence.reset_index(drop=True)
	return outSequence

def main(fileName):
	print 'Running ElementLibrary'
	btFeatures = get_bedtools_features(fileName)
	subsetFeatures = collect_coordinates_for_element_positions(btFeatures,fileName)
	rangeFeatures = get_fasta_for_element_coordinates(subsetFeatures)
	percentage_at_for_element(rangeFeatures['feature'],fileName)
	return rangeFeatures

if __name__ == "__main__":
	main()
