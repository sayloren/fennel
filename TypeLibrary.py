"""
Script to separate by type

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
import pandas as pd
from MethylationLibrary import collect_methylation_data_by_element
from FangsLibrary import run_sliding_window_for_each_nucleotide_string
import GlobalVariables

# do all the analysis for each type
def perform_data_frame_collections_for_group(boolType,fileName):
	if boolType is None or len(boolType) == 0:
		typeWindow, typeNames = None, None
	else:
		typeWindow, typeNames = run_sliding_window_for_each_nucleotide_string(boolType['combineString'],boolType['id'])
	if any(x in GlobalVariables.graphs for x in ['methylation','cluster','methextend']):
		pdMeth = collect_methylation_data_by_element(boolType)
	else: 
		pdMeth = None
	return pdMeth,typeWindow,typeNames

def main(boolType,fileName):
	pdMeth,typeWindow,typeNames = perform_data_frame_collections_for_group(boolType,fileName)
	return pdMeth,typeWindow,typeNames

if __name__ == "__main__":
	main() 
