"""
Script to return tables 

Wren Sayor
July 13 2017

To Do:
Need to calculate the total cpn across the entire genome!!!
Put each tissue type on separate page
Deprecated, at the moment
"""

import argparse
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Make Methylation graphs
def graphTable(TableData,ranTableData,fileName):
	sns.set_style('ticks')
	gs = gridspec.GridSpec(3,1,height_ratios=[3,1,1])
	gs.update(hspace=.8)
	info = str(fileName) # + ', '+ str(len(dfWindow)) + ' - ' "UCES"
	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Table_{0}.pdf'.format(fileName))

	gs = gridspec.GridSpec(1,1,height_ratios=[1])
	gs.update(hspace=.5)

	ax0 = plt.subplot(gs[0,:])
	ax0.set_frame_on(False)
	ax0.set_yticks([])
	ax0.set_xticks([])
# 	ylabels0 = outTable1.columns.str.replace('.bed_PosMethContext','')
	Table0 = ax0.table(cellText=TableData.values,rowLabels=TableData.index,cellLoc='center',rowLoc='center',loc='center')#,colLabels=ylabels1
	ax0.set_title('Plus Strand Methylation',size=8)
	Table0.set_fontsize(8)
	
	sns.despine()
	pp.savefig()
	
	gs = gridspec.GridSpec(1,1,height_ratios=[1])
	gs.update(hspace=.5)
	
	ax1 = plt.subplot(gs[1,:])
	ax1.set_frame_on(False)
	ax1.set_yticks([])
	ax1.set_xticks([])
# 	ylabels1 = outTable2.columns.str.replace('.bed_NegMethContext','')
	Table1 = ax1.table(cellText=ranTableData.values,rowLabels=ranTableData.index,cellLoc='center',rowLoc='center',loc='center')#,colLabels=ylabels2
	ax1.set_title('Minus Strand Methylation',size=8)
	Table1.set_fontsize(8)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	pp.close()

def main(TableData,ranTableData,fileName):
	graphTable(TableData,ranTableData,fileName)

if __name__ == "__main__":
	main()