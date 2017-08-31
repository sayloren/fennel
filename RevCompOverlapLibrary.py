"""
Script to do the reverse complement sorting for the aligned overlap regions

Wren Saylor
August 10 2017

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
from RevCompLibrary import sort_methylation_by_directionality
from RevCompLibrary import sort_sliding_window_by_directionality
import GlobalVariables

# RCsort and return methylation and sliding window computations
def run_directionality_for_overlaps(posStr,negStr):
	compWindow, compNames = sort_sliding_window_by_directionality(negStr,posStr)
	if any(x in GlobalVariables.graphs for x in ['methylation','cluster','methextend']):
		groupMeth = sort_methylation_by_directionality(negStr,posStr)
	else: 
		groupMeth = None
	return groupMeth,compWindow,compNames

def main(posStr,negStr):
	groupMeth,compWindow,compNames = run_directionality_for_overlaps(posStr,negStr)
	return groupMeth,compWindow,compNames

if __name__ == "__main__":
	main()
