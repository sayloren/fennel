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
			t = y+1.0 # include 0
			k = t**2.0 # probability space
			equal.append(t/k)
		else:
			equal.append(0)
	
	sns.set_style('ticks')
	plt.plot(yrange,equal,linewidth=2,color='#3e1638',alpha=0.9)
	plt.axvline(x=n,linewidth=.05,linestyle=':',color='#e7298a')
	plt.axvspan(n-1,n+1,facecolor='#e7298a',label='{0} Bins'.format(n),alpha=0.1)
	plt.xlabel('Bin Size',size=12)
	plt.ylabel('Probability',size=12)
	plt.title('Probability is {0} for Getting Equal Boundary Values for {1} Bins'.format(round(equal[n],2),n),size=18)
	# "{:.2E}".format(Decimal('40800000000.00000000000000'))
	sns.despine()
	plt.savefig('EqualBinProbability.pdf')

def main(n):
	graphComb(n)

if __name__ == "__main__":
	main()