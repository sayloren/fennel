"""
Script to print all the random regions to the same plot, along with element,
run separetly from pipeline, but uses pipeline functions

Wren Saylor
August 18 2017
"""
import argparse
import ElementLibrary
import DirectionLibrary
import FangsLibrary
import MethylationLibrary
import RevCompLibrary
import TypeLibrary
import OverlapLibrary
import BinLibrary
import GraphFangLibrary
import GraphMethLibrary
import GraphSignalLibrary
import BokehLibrary
import GraphTableLibrary
import GraphClusterLibrary
import GraphDendrogramLibrary
import GraphKMeansLibrary
import pandas as pd
import RevCompOverlapLibrary
import BinLibrary
import GlobalVariables
from GraphFangLibrary import collect_sum_two_nucleotides
from FangsLibrary import run_sliding_window_for_each_nucleotide_string
from RevCompLibrary import sort_sliding_window_by_directionality
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev
import scipy.stats as ss
from scipy.stats import mstats
import seaborn as sns

# Arguments
def get_args():
	# File lists
	parser = argparse.ArgumentParser(description="Description")
	parser.add_argument("efile", type=argparse.FileType('rU'), help='A file containing a list of paths to the element files with unique names separated by newlines')
	parser.add_argument("mfile", type=argparse.FileType('rU'), help="A file containing a list of paths to the methylation files with unique names separated by newlines, with data for methylation position (chr, start,stop) and methylation % as fourth column'")
	parser.add_argument("rfile",type=argparse.FileType('rU'), help="A file containing a list of paths to the random regions equable with your elements to plot in contrast")

	# Genome Files
	parser.add_argument("-g", "--genome", type=str, default="hg19.genome")
# 	parser.add_argument("-n", "--nucleosome", type=str, help="A bedgraph file with data for nucleosome positions, form 'chr, start, stop, occupancy'")
# 	parser.add_argument("-s", "--snp", type=str, help="A file with data for snps, form 'chr, start, stop(start+size alt-mutated), ref, ref_size, alt, alt_size, af_adj'")
	parser.add_argument("-fa", "--fasta", type=str, default="hg19.fa")
	parser.add_argument("-o", "--overlapingelements", type=str, default="hg19_0based_exons.bed",help='A file containing all elements to check for overlaps in format "chr start stop", using for exons, but can be anything')#direction

	# Integer Parameters
	parser.add_argument("-t", "--total", type=int, default="600", help='total size of region to look at (region + flanks), should be an even number, suggested to be at least triple your element')
	parser.add_argument("-e", "--element", type=int, default="200", help='size of your element (region without flanks), should be an even number')
	parser.add_argument("-i", "--inset", type=int, default="50", help='size into your element from the boundaries, should be an even number')
	parser.add_argument("-w", "--window", type=int, default="11", help='size of sliding window, should be an odd number, previous studies have used 11')
	parser.add_argument("-b", "--bin", type=int, default="30", help='size of bins used to compare element ends and determine directionality')
	parser.add_argument("-mc", "--thresholdcoverage", type=int, default="10", help='size to threshold uncapped coverage of methylation data to send to % methylation, field often uses 10')
	parser.add_argument("-mp", "--thresholdpercentage", type=int, default="0", help='size to threshold % methylation data')
	parser.add_argument("-mf", "--methylationflank", type=int, default="20", help='The number of base pairs to look at outside of the element for the methylation clusterplots')

	# Specify which groups and graphs to run
	parser.add_argument('-type', "--elementype", default=['all'], nargs='*', choices=['all','intronic','exonic','intergenic'],help='which group types of element to run')
	parser.add_argument('-dir', "--elementdirection", default=[], nargs='*', choices=['+','-','='], help='which group direction of element to run')
	parser.add_argument('-rc', "--reversecomplement",action='store_true', help='if reverse complement sorting required')
	parser.add_argument('-p',"--plots",default=[],nargs='*',choices=['fang','methylation','signal','interactive','cluster','dendrogram','kmean'],help='the available graphs to plot')
	parser.add_argument('-nuc',"--nucleotideline",default=['A','T'],nargs='+',help='type the nucleotide string combinations to search for in the element')
	parser.add_argument('-str',"--stringname",type=str,help='string to add to the outfile name')
	parser.add_argument('-align', "--elementalign",action='store_true', help='if want to align by exonic/intronic crossover')

	# Add additional descriptive file name information
	return parser.parse_args()

# Line graphs
def plot_line_graphs(elementDF,elementRC,collectDF,collectRC,allNames,fileName,collectFile):
	# Get group, mean and standard deviation for AT
	
	# Plot settings
	sns.set_style('ticks')
	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.8) # setting the space between the graphs
	pp = PdfPages('Fangs_{0}.pdf'.format(fileName))
	plt.figure(figsize=(7,7))

	sns.set_palette("husl",n_colors=8)#(len(nucLine)*2)
	ATgroupElement,ATmeanElement,ATstdElement = collect_sum_two_nucleotides(elementDF,allNames,'A','T')
	ATgroupElementRC,ATmeanElementRC,ATstdElementRC = collect_sum_two_nucleotides(elementRC,allNames,'A','T')

	# Plot the mean AT content with a std of 1
	ax0 = plt.subplot(gs[0])
	ax0.plot(GlobalVariables.fillX,ATmeanElement,linewidth=2,label='UCEs')
	for dfNuc,file in zip(collectDF,collectFile):
		ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfNuc,allNames,'A','T')
		ax0.plot(GlobalVariables.fillX,ATmean,linewidth=1,alpha=0.3)#,label='{0}'.format(file)
	ax0.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.set_ylabel('% AT Content',size=8)
	ax0.set_xlabel('Position',size=6)
	ax0.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax0.set_title('Mean AT Content',size=8)
	ax0.set_yticks(ax0.get_yticks()[::2])
	plt.xlim(0,GlobalVariables.num)

	# Plot the std = 1
	ax1 = plt.subplot(gs[1],sharex=ax0)
	ax1.plot(GlobalVariables.fillX,ATmeanElementRC,linewidth=2,label='UCEs')
	for dfNuc,file in zip(collectRC,collectFile):
		ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfNuc,allNames,'A','T')
		ax1.plot(GlobalVariables.fillX,ATmean,linewidth=1,alpha=0.3)#,label='{0}'.format(file)
	ax1.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.set_yticks(ax1.get_yticks()[::2])
	ax1.set_xlabel('Position',size=6)
	ax1.set_ylabel('% AT Content',size=8)
	ax1.set_title('Mean AT Content, Reverse Completement Sorted',size=8)
	plt.setp(ax1.get_xticklabels(), visible=True)
	ax1.legend(loc=0,fontsize=5,labelspacing=0.05)
	sns.despine()
	pp.savefig()
	pp.close()

def main():
	# Collect arguments
	args = get_args()
	# Set global variables
	GlobalVariables.main(args)
	collectDF = []
	collectFile = []
	collectRC = []
	
	elementFile = []
	for fileName in GlobalVariables.eFiles:
		elementFile.append(fileName)
		rangeFeatures = ElementLibrary.main(fileName)
		directionFeatures,directionBins = DirectionLibrary.main(rangeFeatures,fileName)
		elementDF, elementNames = run_sliding_window_for_each_nucleotide_string(rangeFeatures['combineString'],rangeFeatures['id'])
		negStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '-')])
		posStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '+')])
		elementRC, elementRCNames = sort_sliding_window_by_directionality(negStr,posStr)

	# for each element file provided
	for fileName in GlobalVariables.rFiles:
		rangeFeatures = ElementLibrary.main(fileName)
		directionFeatures,directionBins = DirectionLibrary.main(rangeFeatures,fileName)
		# separate by direction
		typeWindow, typeNames = run_sliding_window_for_each_nucleotide_string(rangeFeatures['combineString'],rangeFeatures['id'])
		collectDF.append(typeWindow)
		collectFile.append(fileName)
		negStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '-')])
		posStr = (directionFeatures[(directionFeatures['compareBoundaries'] == '+')])
		compWindow, compNames = sort_sliding_window_by_directionality(negStr,posStr)
		collectRC.append(compWindow)
	paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}_{6}'.format(GlobalVariables.uce,GlobalVariables.inuce,GlobalVariables.num,GlobalVariables.binDir,GlobalVariables.window,elementFile,GlobalVariables.stringName)
	plot_line_graphs(elementDF,elementRC,collectDF,collectRC,typeNames,'all_random_{0}'.format(paramlabels),collectFile)

# 			# By Type
# 			for type in GlobalVariables.typeList:
# 				typeBool,typeMeth,typeWindow,typeNames = separate_dataframe_by_group(type,directionFeatures,'type',fileName)
# 				if GlobalVariables.revCom:
# 					typercMeth,typercWindow,typercNames = RevCompLibrary.main(typeBool)

if __name__ == "__main__":
	main()

"""
                                 Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "{}"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright {yyyy} {name of copyright owner}

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