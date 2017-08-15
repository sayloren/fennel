"""

Script to determine how the probability of getting boundaries with the same
value changes as the bin size of the element boundaries changes

Wren Saylor
June 20 2017

"""

import argparse
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb
import seaborn as sns

def runProbability(yrange,prange):
	equal = []
	for y,p in zip(yrange,prange):
		if y >= 0:
			t = y + 1.0 # include 0% don't think this is needed
			k = t**2.0 # probability space
			equal.append(t/k)
		else:
			equal.append(0)
	return equal

def graphComb(n,emp,ranemp,paramlabels):
# https://betterexplained.com/articles/understanding-the-birthday-paradox/
	yrange = np.arange(1,n*2)
	prange = np.empty(n*2)
	prange.fill(0.5)
	equal = runProbability(yrange,prange)
	
	sns.set_style('ticks')
	sns.set_palette("husl",n_colors=8)#(len(nucLine)*2)

	plt.figure(figsize=(3.5,3.5))
	plt.plot(yrange,equal,linewidth=2,alpha=0.9,label='Theoretical')
	plt.plot(yrange,emp['Equal'],linewidth=2,alpha=0.9,label='Observed')
	plt.plot(yrange,ranemp['Equal'],linewidth=2,alpha=0.9,label='Random')
# 	plt.fill_between(yrange,emp['Min'],emp['Max'],alpha=0.2)#,label='Observed Range'
# 	plt.fill_between(yrange,ranemp['Min'],ranemp['Max'],alpha=0.2)#,label='Random Range'
	plt.axvline(x=n,linewidth=.05,linestyle=':')#label='{:0.1e}'.format(equal[n])
	plt.axvspan(n-1,n+1,alpha=0.1)#,facecolor='#e7298a'
	plt.xlabel('Bin Size',size=12)
	plt.ylabel('Probability',size=12)
	plt.title('Equal Boundary for {0} Bins'.format(n),size=16)
	plt.legend(loc=0,fontsize=6,labelspacing=0.05)
# 	plt.text(n-9,.5,'{:0.1e}'.format(equal[n]),size=12,clip_on=False)
	plt.tight_layout()
	sns.despine()
	plt.savefig('Probability_{0}.pdf'.format(paramlabels))

def main(n,emp,ranemp,paramlabels):
	print 'Running BinLibrary'
	graphComb(n,emp,ranemp,paramlabels)

if __name__ == "__main__":
	main()