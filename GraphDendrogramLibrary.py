"""
Script to graph the dendrogram connections

Wren Saylor
August 4 2017

To Do:
The arm colors don't work yet, need to change thresholding

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
from fastcluster import linkage
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.colors import rgb2hex,colorConverter,ListedColormap
import scipy
from scipy import cluster
from scipy.cluster.hierarchy import dendrogram,set_link_color_palette
from scipy.spatial import distance
from scipy.cluster import hierarchy
from GraphFangLibrary import collect_sum_two_nucleotides
import GlobalVariables

# Make cluster classes based on the den values
def find_cluster_classification(den):
	#https://stackoverflow.com/questions/27924813/extracting-clusters-from-seaborn-clustermap
	#http://www.nxn.se/valent/extract-cluster-elements-by-color-in-python
	label='ivl'
	cluster_idxs = defaultdict(list)
	for c, pi in zip(den['color_list'],den['icoord']):
		for leg in pi[1:3]:
			i = (leg - 5.0)/10.0
			if abs(i-int(i)) < 1e5:
				cluster_idxs[c].append(int(i))
	cluster_classes = {}
	for c, l in cluster_idxs.items():
		i_l = [den[label][i] for i in l]
		cluster_classes[c] = i_l
	return cluster_classes

# Make a column for the cluster class values
def make_cluster_column(ATelement,Elementclusters):
	cluster = []
	for i in ATelement.T.index:
		included=False
		for j in Elementclusters.keys():
			if i in Elementclusters[j]:
				cluster.append(j)
				included=True
		if not included:
			cluster.append(None)
	return cluster

# Make some graphs for fangs
def graph_dendrogram_branches(dfWindow,ranWindow,names,fileName):

	# Get group, mean and standard deviation for AT
	ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(ranWindow,names,'A','T')
	ATelement = ATgroup.T[(GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank):(GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank)]
	ranATelement = ranATgroup.T[(GlobalVariables.plotLineLocationThree-GlobalVariables.methylationflank):(GlobalVariables.plotLineLocationFour+GlobalVariables.methylationflank)]
	print 'Extracted just element and methylation flank, size {0}'.format(len(ATelement))

	# Title info
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Dendrogram_{0}.pdf'.format(fileName))
	sns.set_palette("husl",n_colors=8)
	palette = sns.color_palette()
	set_link_color_palette(map(rgb2hex, palette))
	sns.set_style('white')

	#http://nbviewer.jupyter.org/gist/vals/150ec97a5b7db9c82ee9
	link = linkage(ATelement.T)
	plt.figure(figsize=(100,10))
	den = dendrogram(link,labels=ATelement.T.index,leaf_font_size=3,color_threshold='#AAAAAA')#,link_color_func=lambda x:___[x] Have to make dict of len df, with clustered colors
	plt.xticks(rotation=90,fontsize=8)
	sns.despine()

	plt.tight_layout()

	sns.despine()
	pp.savefig()
	pp.close()

def main(dfWindow,ranWindow,names,fileName):
	print 'Running graph_dendrogram_branchesLibrary'
	graph_dendrogram_branches(dfWindow,ranWindow,names,fileName)

if __name__ == "__main__":
	main()
