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
	sns.set_palette("husl",n_colors=8)#(len(nucLine)*2)

	plt.figure(figsize=(3.5,3.5))
	plt.plot(yrange,equal,linewidth=2,alpha=0.9)#,color='#3e1638'
	plt.axvline(x=n,linewidth=.05,linestyle=':')#,color='#e7298a',label='{:0.1e}'.format(equal[n])
	plt.axvspan(n-1,n+1,alpha=0.1)#,facecolor='#e7298a'
	plt.xlabel('Bin Size',size=12)
	plt.ylabel('Probability',size=12)
	plt.title('Equal Boundary for {0} Bins'.format(n),size=16)
# 	plt.legend(loc=2,fontsize=12,labelspacing=0.05)
	plt.text(n-9,.9,'{:0.1e}'.format(equal[n]),size=12,clip_on=False)
	plt.tight_layout()
	sns.despine()
	plt.savefig('Probability.pdf')

def main(n):
	graphComb(n)

if __name__ == "__main__":
	main()