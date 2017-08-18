"""
Script to call and execute different sets of analyses for Element Structure

Wren Saylor
July 5 2017

To Do:
k means / scikit learn
dendrogram color clusters (what clustering parameters are used in literature, how do uces compare to other clustered groups)
cluster by difference between boundaries
snps, nucleosomes, replication origins/forks, tss

AT balance random
exons - split intron/intergenic, include exon direcitonality, move cross boundary to tss, rc sorting
which are consitently = 

unittest
Apache License 2.0

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

# set command line arguments
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

# for type and direction, separate the groups and run the analyses
def separate_dataframe_by_group(List,directionFeatures,typecolumn,fileName):
	# subset by bool presence
	bool = (directionFeatures[directionFeatures[typecolumn] == List])
	# if there is nothing in that set, skip
	if len(bool.index) != 0:
		Meth,dfWindow,Names = TypeLibrary.main(bool,fileName)
	return bool,Meth,dfWindow,Names

# the plotting options, if in the list of plot flags, run graph
def plot_graph_by_arg(pdMeth,rnMeth,dfWindow,names,ranWindow,fileName):
	if 'fang' in GlobalVariables.graphs:
		GraphFangLibrary.main(dfWindow,names,ranWindow,fileName)
	if 'signal' in GlobalVariables.graphs:
		GraphSignalLibrary.main(dfWindow,names,ranWindow,fileName)
	if 'methylation' in GlobalVariables.graphs:
		GraphMethLibrary.main(pdMeth,rnMeth,fileName)
	if 'interactive' in GlobalVariables.graphs:
		BokehLibrary.main(dfWindow,ranWindow,fileName)
	if 'cluster' in GlobalVariables.graphs:
		GraphClusterLibrary.main(dfWindow,ranWindow,pdMeth,rnMeth,names,fileName)
	if 'dendrogram' in GlobalVariables.graphs:
		GraphDendrogramLibrary.main(dfWindow,ranWindow,names,fileName)
	if 'kmean' in GlobalVariables.graphs:
		GraphKMeansLibrary.main(dfWindow,ranWindow,names,fileName)

# make a table for element and random region
def plot_chi_square_table(TableData,Title,ranTableData,ranTitle,fileName):
	GraphTableLibrary.main(TableData,Title,ranTableData,ranTitle,fileName)

# collect the counts for how many of each boundary direction for each type
def collect_counts_for_element_type(rangeFeatures):
	grouptype = rangeFeatures.groupby('type')['compareBoundaries'].value_counts()
	pdgroup = pd.DataFrame(grouptype)
	pdgroup.columns = ['counts']
	pvtable = pd.pivot_table(pdgroup,index=['type'],columns=['compareBoundaries'],values=['counts'])
	pvtable.columns.name = None
	pvtable.index.name = None
	pvtable.loc['all']= pvtable.sum()
	pvout = pvtable.xs('counts', axis=1, drop_level=True)
	return pvout

def main():
	# Collect arguments
	args = get_args()
	
	# Set global variables
	GlobalVariables.main(args)

	# for each element file provided
	for fileName in GlobalVariables.eFiles:
		# get the region info to work with
		rangeFeatures = ElementLibrary.main(fileName)
		# separate by direction
		directionFeatures,directionBins = DirectionLibrary.main(rangeFeatures,fileName)
		
		# for each random file provided
		for randomFile in GlobalVariables.rFiles:
			# Collect out filename labels 
			paramlabels = '{0}_{1}_{2}_{3}_{4}_{5}_{6}_{7}'.format(GlobalVariables.uce,GlobalVariables.inuce,GlobalVariables.num,GlobalVariables.binDir,GlobalVariables.window,fileName,randomFile,GlobalVariables.stringName)
			
			# get the region info to work with
			randomFeatures = ElementLibrary.main(randomFile)
			# separate by direction
			randirFeatures,randirBins = DirectionLibrary.main(randomFeatures,randomFile)
			
			# Plot boundary probabilities
			BinLibrary.main(directionBins,randirBins,paramlabels)

			# Make table for the count of each direction for each type element
			elementGroup = collect_counts_for_element_type(rangeFeatures)
			randomGroup = collect_counts_for_element_type(randirFeatures)
			plot_chi_square_table(elementGroup,'Elements',randomGroup,'Random Regions','Direction_Count_{0}'.format(paramlabels))
			
			# All elements
			if 'all' in GlobalVariables.typeList:
				GlobalVariables.typeList.remove('all')
				pdMeth,allWindow,allNames = TypeLibrary.main(rangeFeatures,fileName)
				rnMeth,ranWindow,ranNames = TypeLibrary.main(randomFeatures,fileName)
				plot_graph_by_arg(pdMeth,rnMeth,allWindow,allNames,ranWindow,'all_{0}'.format(paramlabels))
				if GlobalVariables.revCom:
					revMeth,revWindow,revNames = RevCompLibrary.main(directionFeatures)
					ranrevMeth,ranrevWindow,ranrevNames = RevCompLibrary.main(randirFeatures)
					plot_graph_by_arg(revMeth,ranrevMeth,revWindow,revNames,ranrevWindow,'revComp_all_{0}'.format(paramlabels))

# 			# By Type
			for type in GlobalVariables.typeList:
				typeBool,typeMeth,typeWindow,typeNames = separate_dataframe_by_group(type,directionFeatures,'type',fileName)
				rantypeBool,rantypeMeth,rantypeWindow,rantypeNames = separate_dataframe_by_group(type,randirFeatures,'type',randomFile)
				plot_graph_by_arg(typeMeth,rantypeMeth,typeWindow,typeNames,rantypeWindow,'{0}_{1}'.format(type,paramlabels))
				if GlobalVariables.revCom:
					typercMeth,typercWindow,typercNames = RevCompLibrary.main(typeBool)
					rantypercMeth,rantypercWindow,rantypercNames = RevCompLibrary.main(rantypeBool)
					plot_graph_by_arg(typercMeth,rantypercMeth,typercWindow,typercNames,rantypercWindow,'revComp_{0}_{1}'.format(type,paramlabels))

			# By Direction
			for dir in GlobalVariables.dirList:
				dirBool,dirMeth,dirWindow,dirNames = separate_dataframe_by_group(dir,directionFeatures,'compareBoundaries',fileName)
				randirBool,randirMeth,randirWindow,randirNames = separate_dataframe_by_group(dir,randirFeatures,'compareBoundaries',randomFile)
				plot_graph_by_arg(dirMeth,randirMeth,dirWindow,dirNames,randirWindow,'all_{0}_{1}'.format(dir,paramlabels))
				for type in GlobalVariables.typeList:
					typeBool,typeMeth,typeWindow,typeNames = separate_dataframe_by_group(type,directionFeatures,'type',fileName)
					rantypeBool,rantypeMeth,rantypeWindow,rantypeNames = separate_dataframe_by_group(type,randirFeatures,'type',randomFile)
					plot_graph_by_arg(typeMeth,rantypeMeth,typeWindow,typeNames,rantypeWindow,'{0}_{1}_{2}'.format(type,dir,paramlabels))

			# Re-align by exon/intron crossover
			if GlobalVariables.overlapInset:
				ranstartcrossboundary,ranendcrossboundary,rancompleteelement,raninteriorelement,ranoverlapTable = OverlapLibrary.main(randomFeatures,fileName)
				rancrossMeth,rancrossWindow,rancrossNames  = RevCompOverlapLibrary.main(ranstartcrossboundary,ranendcrossboundary)
				startcrossboundary,endcrossboundary,completeelement,interiorelement,overlapTable = OverlapLibrary.main(rangeFeatures,fileName)
				crossMeth,crossWindow,crossNames  = RevCompOverlapLibrary.main(startcrossboundary,endcrossboundary)
				completeelement.name = 'WithinExon'
				interiorelement.name = 'ContainsExon'
				plot_chi_square_table(overlapTable.T,'Exonic Elements',ranoverlapTable.T,'Exonic Random Regions','Exonic_overlap_Count_{0}'.format(paramlabels))
				plot_graph_by_arg(crossMeth,rancrossMeth,crossWindow,crossNames,rancrossWindow,'align_Crossboundary_{0}'.format(paramlabels))
				for element,random in zip([completeelement,interiorelement],[rancompleteelement,raninteriorelement]):
					alignMeth,alignWindow,alignNames = TypeLibrary.main(element,fileName)
					ranalignMeth,ranalignWindow,ranalignNames = TypeLibrary.main(random,fileName)
					plot_graph_by_arg(alignMeth,ranalignMeth,alignWindow,alignNames,ranalignWindow,'align_{0}_{1}'.format(element.name,paramlabels))

# 					if GlobalVariables.revCom: # maybe use group separate?
# 						aligndirFeatures = DirectionLibrary.main(element,fileName)
# 						alignranFeatures = DirectionLibrary.main(random,fileName)
# 						alignrcMeth,alignrcWindow,alignrcNames = RevCompLibrary.main(aligndirFeatures)
# 						ranalignrcMeth,ranalignrcWindow,ranalignrcNames = RevCompLibrary.main(alignranFeatures)
# 						plot_graph_by_arg(alignrcMeth,ranalignrcMeth,alignrcWindow,alignrcNames,ranalignrcWindow,'align_{0}_{1}'.format(element.name,paramlabels))

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