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
from scipy.stats import chisquare

# Make Methylation graphs
def graphTable(TableData,Title,ranTableData,ranTitle,fileName):
	sns.set_style('ticks')
	info = str(fileName)
# 	plt.suptitle(info,fontsize=10)
	pp = PdfPages('Table_{0}.pdf'.format(fileName))
	gs = gridspec.GridSpec(3,1,height_ratios=[1,1,1])
	gs.update(hspace=.1)
	plt.tight_layout()
	plt.figure(figsize=(12,6))

	ax0 = plt.subplot(gs[0])
	ax0.set_frame_on(False)
	ax0.set_yticks([])
	ax0.set_xticks([])
	Table0 = ax0.table(cellText=TableData.values,rowLabels=TableData.index,colLabels=TableData.columns,cellLoc='center',rowLoc='center',loc='center')
	ax0.set_title(Title,size=14)
	Table0.set_fontsize(14)
	
	ax1 = plt.subplot(gs[1])
	ax1.set_frame_on(False)
	ax1.set_yticks([])
	ax1.set_xticks([])
	Table1 = ax1.table(cellText=ranTableData.values,rowLabels=ranTableData.index,colLabels=ranTableData.columns,cellLoc='center',rowLoc='center',loc='center')
	ax1.set_title(ranTitle,size=14)
	Table1.set_fontsize(14)

	outchisquare = []
	for x,y in zip(TableData,ranTableData):
		outchisquare.append(chisquare(TableData[x],ranTableData[y]))
	pdchi = pd.DataFrame(outchisquare)
	pdchi.index = TableData.T.index
	ax2 = plt.subplot(gs[2])
	ax2.set_frame_on(False)
	ax2.set_yticks([])
	ax2.set_xticks([])
	Table2 = ax2.table(cellText=pdchi.values.round(3),rowLabels=pdchi.index,colLabels=pdchi.columns,cellLoc='center',rowLoc='center',loc='center')
	ax2.set_title('Chi Square Test',size=14)
	Table2.set_fontsize(14)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	pp.close()

def main(TableData,Title,ranTableData,ranTitle,fileName):
	graphTable(TableData,Title,ranTableData,ranTitle,fileName)

if __name__ == "__main__":
	main()