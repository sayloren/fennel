"""

Script to determine how the probability of getting boundaries with the same
value changes as the bin size of the element boundaries changes

Wren Saylor
June 20 2017

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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy
import seaborn as sns
import math
import GlobalVariables

# Compute the theoretical probability
def run_boundary_probability_calculation(yrange,p):
	# https://math.stackexchange.com/questions/151810/probability-of-3-heads-in-10-coin-flips
	# need to include probability different than 0.5
	equal = []
	for y in yrange:
		totalPerm = pow(2,y)
		if y >= 0:
			permsum = []
			for k in range (0, (y + 1)):
				permuationK = math.factorial(y)/(math.factorial(k)*(math.factorial(y -k)))#*pow(p,k)*pow((1-p),y-k)
				floatPForKHeads = float(permuationK)/float(totalPerm)
				floatPForKHeadsBothSides = floatPForKHeads * floatPForKHeads
				permsum.append(floatPForKHeadsBothSides)
			probabilitiesfloatPEqualSides = sum(permsum)
			equal.append(probabilitiesfloatPEqualSides)
		else:
			equal.append(0)
	return equal

# Graph the probability
def graph_equal_boundary_probability(emp,ranemp,paramlabels):
	yrange = numpy.arange(1,GlobalVariables.binDir*2)
	equal = run_boundary_probability_calculation(yrange,0.5)
	
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
	plt.text(GlobalVariables.binDir-9,.5,'{:0.1e}'.format(equal[GlobalVariables.binDir]),size=12,clip_on=False)
	plt.tight_layout()
	sns.despine()
	plt.savefig('Probability_{0}.pdf'.format(paramlabels))

def main(emp,ranemp,paramlabels):
	print 'Running BinLibrary'
	graph_equal_boundary_probability(emp,ranemp,paramlabels)

if __name__ == "__main__":
	main()
