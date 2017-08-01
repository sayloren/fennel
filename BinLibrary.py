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

def graphComb(n):
	# called birthday problem
	yrange = np.arange(0,n*2)
	equal = []
	for y in yrange:
		if y >= 0:
			t = y+1.0 # include 0%
			k = t**2.0 # probability space
			equal.append(t/k)
		else:
			equal.append(0)
	
	sns.set_style('ticks')
	plt.figure(figsize=(8,8))
	plt.plot(yrange,equal,linewidth=2,color='#3e1638',alpha=0.9,label='Probability is {:0.1e}'.format(equal[n]))
	plt.axvline(x=n,linewidth=.05,linestyle=':',color='#e7298a')
	plt.axvspan(n-1,n+1,facecolor='#e7298a',label='{0} Bins'.format(n),alpha=0.1)
	plt.xlabel('Bin Size',size=12)
	plt.ylabel('Probability',size=12)
	plt.title('Equal Boundary for {0} Bins'.format(n),size=18)
	plt.legend(loc=0,fontsize=12,labelspacing=0.05)
	sns.despine()
	plt.savefig('Probability.pdf')

def main(n):
	graphComb(n)

if __name__ == "__main__":
	main()