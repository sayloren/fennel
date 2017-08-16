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
import math
import GlobalVariables

def runProbability(yrange):
	# https://math.stackexchange.com/questions/151810/probability-of-3-heads-in-10-coin-flips
	# need to include probability different than 0.5
	equal = []
	for y in yrange:
		totalPerm = pow(2,y)
		if y >= 0:
			permsum = []
			for k in range (0, (y + 1)):
				permuationK = math.factorial(y)/(math.factorial(k)*(math.factorial(y -k)))
				floatPForKHeads = float(permuationK)/float(totalPerm)
				floatPForKHeadsBothSides = floatPForKHeads * floatPForKHeads
				permsum.append(floatPForKHeadsBothSides)
			probabilitiesfloatPEqualSides = sum(permsum)
			equal.append(probabilitiesfloatPEqualSides)
		else:
			equal.append(0)
	return equal

def graphComb(emp,ranemp,paramlabels):
	yrange = np.arange(1,GlobalVariables.binDir*2)
	equal = runProbability(yrange)
	
	sns.set_style('ticks')
	sns.set_palette("husl",n_colors=8)

	plt.figure(figsize=(3.5,3.5))
	plt.plot(yrange,equal,linewidth=2,alpha=0.9,label='Theoretical')
	plt.plot(yrange,emp['Equal'],linewidth=2,alpha=0.9,label='Observed')
	plt.plot(yrange,ranemp['Equal'],linewidth=2,alpha=0.9,label='Random')
	plt.axvline(x=GlobalVariables.binDir,linewidth=.05,linestyle=':')#label='{:0.1e}'.format(equal[n])
	plt.axvspan(GlobalVariables.binDir-1,GlobalVariables.binDir+1,alpha=0.1)
	plt.xlabel('Bin Size',size=12)
	plt.ylabel('Probability',size=12)
	plt.title('Equal Boundary for {0} Bins'.format(GlobalVariables.binDir),size=16)
	plt.legend(loc=0,fontsize=6,labelspacing=0.05)
# 	plt.text(n-9,.5,'{:0.1e}'.format(equal[n]),size=12,clip_on=False)
	plt.tight_layout()
	sns.despine()
	plt.savefig('Probability_{0}.pdf'.format(paramlabels))

def main(emp,ranemp,paramlabels):
	print 'Running BinLibrary'
	graphComb(emp,ranemp,paramlabels)

if __name__ == "__main__":
	main()