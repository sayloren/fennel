"""

Script to run the same analysises as Fang element library, but makes bokeh plots

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
from bokeh.models import BoxSelectTool, LassoSelectTool, Spacer
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, curdoc
from bokeh.plotting import figure, output_file, show
from bokeh.layouts import row, column
from bokeh.layouts import widgetbox
from bokeh.models.widgets import Select
from bokeh.models.widgets import Toggle
import pandas as pd
import seaborn
from GraphFangLibrary import collect_sum_two_nucleotides
import GlobalVariables

# Make interactive plots
def graph_interactive_bokeh(dfWindow,ranWindow,fileName):
	seaborn.set_palette("husl",n_colors=8)

	# Get group, mean and standard deviation for AT
	ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(ranWindow,names,'A','T')

	source = ColumnDataSource(data=dict(x=GlobalVariables.fillX,mean=ATgroup.mean(),std=ATgroup.std(),rmean=ranATgroup.mean(),rstd=ranATgroup.std()))
	output_file('Interactive_{0}.html'.format(fileName))
	p = figure(plot_width=1500, plot_height=600, min_border=10, min_border_left=50,toolbar_location="above",title="Mean AT Content Across Base Pair Position")
	p.line('x','mean',line_width=2,source=source)
	p.yaxis.axis_label = "% AT Content"
	p.xaxis.axis_label = "Nucleotide Postion"
	p.background_fill_color = "#fafafa"
	
	select = Select(title="Option:", value="All", option=["All", "Exonic", "Intronic", "Intergenic"])
	#toggle = Toggle(label="Reverse Complement")
	
	sd = figure(plot_width=1500, plot_height=200, x_range=p.x_range, min_border=10, min_border_left=50,title="Standard Deviation")
	sd.line('x','std',line_width=2,color='#3e1638',alpha=.5,source=source)
	sd.background_fill_color = "#fafafa"
	
	widgets = row(select,toggle)
	show(column(widgets,p,sd))

def main(dfWindow,ranWindow,fileName):
	print 'Running BokehLibrary'
	graph_interactive_bokeh(dfWindow,ranWindow,fileName)

if __name__ == "__main__":
	main()
