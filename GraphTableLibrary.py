"""
Script to return tables 

Wren Sayor
July 13 2017
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
def graphTable(TableData,Title,ranTableData,ranTitle,fileName):
	sns.set_style('ticks')
	info = str(fileName)
# 	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Table_{0}.pdf'.format(fileName))
	plt.tight_layout()
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.5)

	ax0 = plt.subplot(gs[0])
	ax0.set_frame_on(False)
	ax0.set_yticks([])
	ax0.set_xticks([])
	Table0 = ax0.table(cellText=TableData.values,rowLabels=TableData.index,colLabels=TableData.columns,cellLoc='center',rowLoc='center',loc='center')
	ax0.set_title(Title,size=8)
	Table0.set_fontsize(8)
	
	ax1 = plt.subplot(gs[1])
	ax1.set_frame_on(False)
	ax1.set_yticks([])
	ax1.set_xticks([])
	Table1 = ax1.table(cellText=ranTableData.values,rowLabels=ranTableData.index,colLabels=ranTableData.columns,cellLoc='center',rowLoc='center',loc='center')
	ax1.set_title(ranTitle,size=8)
	Table1.set_fontsize(8)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	pp.close()

def main(TableData,Title,ranTableData,ranTitle,fileName):
	graphTable(TableData,Title,ranTableData,ranTitle,fileName)

if __name__ == "__main__":
	main()