"""
Script to graph Signal analyses

Wren Saylor
July 5 2017

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
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from numpy import sin, linspace, pi
import numdifftools as nd
import numpy as np
import pandas as pd
import scipy.stats as ss
import scipy as sp
import scipy.fftpack
from scipy.interpolate import splrep, splev
from scipy import signal
from scipy.stats import mstats
import seaborn as sns
from GraphFangLibrary import collect_sum_two_nucleotides
import GlobalVariables

# Get just the elemenet
def get_element_coordinates(region):
	element = region[GlobalVariables.plotLineLocationThree:GlobalVariables.plotLineLocationFour]
	return element

# Get just the downstream flank
def get_down_stream_flank_coordinates(region):
	dFlank = region[GlobalVariables.plotLineLocationFour:(GlobalVariables.num-GlobalVariables.window-GlobalVariables.window)]
	return dFlank

# Get just the upstream flank
def get_up_stream_flank_coordinates(region):
	uFlank = region[GlobalVariables.window:GlobalVariables.plotLineLocationThree]
	return uFlank

# Perform fourier transform
def perform_fourier_transfroms(region):
	Fs = 1.0 # sampling rate
	Ts = 1.0/Fs # sampling interval
	
	# length of the signal
	nsd = len(region)
	ksd = np.arange(nsd)
	Tsd = nsd/Fs
	
	# two sides frequency range
	frqsd = ksd/Tsd
	
	# one side frequncy range
	frqsd = frqsd[range(nsd/2)]
	
	# fft computing and normalization
	Ysd = np.fft.fft(region)/nsd
	Ysd = Ysd[range(nsd/2)]
	
	return frqsd, Ysd

# Collect each UCEs second derivative, find their inflection points
def behaviorlocate_second_derivative_inflection_pointsUCE(ATgroup):
	inflectionPts = []
	for index, row in ATgroup.iterrows():
		collectUCE = []
		collectUCE.append(index)
		f = splrep(GlobalVariables.fillX,row,k=5,s=11)
		smoothMean = splev(GlobalVariables.fillX,f)
		secondDer = splev(GlobalVariables.fillX,f,der=2)
		secondDer[0:GlobalVariables.window] = 0 # small edge effect
		secondDer[-GlobalVariables.window:] = 0 # small edge effect
		peaks = signal.find_peaks_cwt(secondDer,np.arange(1,45))
		peaksOut = [s for s in peaks if s > (GlobalVariables.window*3) and s < (GlobalVariables.num-(GlobalVariables.window*3))] # Get rid of edge effects
		collectUCE.append(peaksOut)
		inflectionPts.append(collectUCE)
	peaksUCE = pd.DataFrame(inflectionPts)
	print 'Collected inflection points for {0} elements'.format(len(ATgroup.index))
	return peaksUCE

# Get smoothed mean, first and second derivatives
def collect_smoothed_lines(ATmean):
	f = splrep(GlobalVariables.fillX,ATmean,k=5,s=11)
	smoothMean = splev(GlobalVariables.fillX,f)
	firstDer = splev(GlobalVariables.fillX,f,der=1)
	firstDer[0:GlobalVariables.halfwindow] = 0 # small edge effect
	firstDer[-GlobalVariables.halfwindow:] = 0 # small edge effect
	secondDer = splev(GlobalVariables.fillX,f,der=2)
	secondDer[0:GlobalVariables.window] = 0 # small edge effect
	secondDer[-GlobalVariables.window:] = 0 # small edge effect
	return smoothMean,firstDer,secondDer

# Get inflection points for seconder derivative
def locate_second_derivative_inflection_points(ATgroup):
	infUCEpeaks = behaviorlocate_second_derivative_inflection_pointsUCE(ATgroup)
	infUCEpeaks.columns = ['id','locate_second_derivative_inflection_points']
	inflectionList = infUCEpeaks['locate_second_derivative_inflection_points'].apply(pd.Series).stack().tolist()
	print 'Collected all inflections points to list for {0} elements'.format(len(ATgroup.index))
	return inflectionList, infUCEpeaks

# Make signal graphs
def graph_signal_lines(dfWindow,names,ranWindow,fileName):
	
	# Get group, mean and standard deviation for AT
	ATgroup,ATmean,ATstd = collect_sum_two_nucleotides(dfWindow,names,'A','T')
	ranATgroup,ranATmean,ranATstd = collect_sum_two_nucleotides(ranWindow,names,'A','T')

	# File name
	info = str(fileName) + ', '+ str(len(ATgroup.index)) + ' - ' "UCES"

	# Plot settings
	sns.set_style('ticks')
	plt.suptitle(info,fontsize=10)
	sns.set_palette("husl",n_colors=8)#(len(nucLine)*2)
	plt.figure(figsize=(7,7))

	# Filename
	pp = PdfPages('Signal_{0}.pdf'.format(fileName))

	# Get smoothed mean, first and second derivatives
	smoothMean,firstDer,secondDer = collect_smoothed_lines(ATmean)
	ransmoothMean,ranfirstDer,ransecondDer = collect_smoothed_lines(ranATmean)

	gs = gridspec.GridSpec(3,3,height_ratios=[1,1,1])
	gs.update(hspace=.65)
	# Plot smoothed mean AT
	ax0 = plt.subplot(gs[0,:])
	ax0.plot(GlobalVariables.fillX,smoothMean,linewidth=1,alpha=0.9,label='Element')
	ax0.plot(GlobalVariables.fillX,ransmoothMean,linewidth=1,alpha=0.9,label='Random')
	ax0.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax0.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax0.set_yticks(ax0.get_yticks()[::2])
	ax0.set_ylabel('% AT Content',size=8)
	ax0.set_title('Fitted Mean AT Content',size=8)
	
	# First derivative
	ax1 = plt.subplot(gs[1,:],sharex=ax0)
	ax1.plot(GlobalVariables.fillX,firstDer,linewidth=1,alpha=0.8,label='Element')
	ax1.plot(GlobalVariables.fillX,ranfirstDer,linewidth=1,alpha=0.8,label='Random')
	ax1.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax1.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax1.axhline(y=0,linewidth=.1,alpha=0.3)#,color='#bd4973'
	ax1.set_yticks(ax1.get_yticks()[::2])
	ax1.set_ylabel('Amplitude',size=8)
	ax1.set_title('First Derivative of Fitted Mean',size=8)
	
	# Second derivative
	ax2 = plt.subplot(gs[2,:],sharex=ax0)
	peakMean = signal.find_peaks_cwt(secondDer,np.arange(1,45)).astype(int)
	print 'Found peaks for elements second derivative'
	ax2.plot(GlobalVariables.fillX,secondDer,linewidth=1,alpha=0.7,label='Element')
	ax2.plot(GlobalVariables.fillX,ransecondDer,linewidth=1,alpha=0.7,label='Random')
	ax2.scatter(peakMean,secondDer[peakMean],marker='.')#,color='#ae3e9e'
	ax2.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax2.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax2.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax2.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax2.axhline(y=0,linewidth=.1,alpha=0.3)#,color='#bd4973'
	ax2.set_ylabel('Amplitude',size=8)
	ax2.set_xlabel('Position',size=6)
	ax2.set_yticks(ax2.get_yticks()[::2])
	ax2.set_title('Second Derivative of Fitted Mean',size=8)
	sns.despine()
	plt.savefig(pp, format='pdf')
	print 'Plotted mean, first and second derivatives for elements and random regions'
	
	# Frequency of inflection points
	gs = gridspec.GridSpec(1,1,height_ratios=[1])
	ax3 = plt.subplot(gs[0])
	inflectionList,infUCEpeaks = locate_second_derivative_inflection_points(ATgroup)
	raninflectionList,raninfUCEpeaks = locate_second_derivative_inflection_points(ranATgroup)
	IFbins = GlobalVariables.num/10
	ax3.hist(inflectionList,IFbins,alpha=0.3,label='Element')
	ax3.hist(raninflectionList,IFbins,alpha=0.3,label='Random')
	ax3.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax3.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax3.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax3.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax3.set_yticks(ax3.get_yticks()[::2])
	ax3.set_ylabel('Frequency',size=8)
	ax3.legend(loc=0,fontsize=5,labelspacing=0.1)
	ax3.set_xlabel('Inflection Point Location',size=8)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	print 'Plotted inflection point frequency for elements and random regions'

	gs = gridspec.GridSpec(2,1,height_ratios=[1,1])
	gs.update(hspace=.65)
	ax5 = plt.subplot(gs[0],sharex=ax0)
	endRange = 25
	widths = np.arange(1, endRange)
	cwtmatr = signal.cwt(firstDer, signal.ricker, widths)
	ax5.imshow(cwtmatr,cmap='RdPu',extent=[0, (GlobalVariables.num-GlobalVariables.window), 1, endRange],aspect='auto',vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())
	ax5.set_xlabel('Position',size=6)
	ax5.set_yticks(ax5.get_yticks()[::2])
	ax5.set_title('Continuous Wavelet Transformation Convolved Over Range {0}-{1} for the First Derivative, Element'.format(widths[0],endRange),size=8)
	
	ax6 = plt.subplot(gs[1],sharex=ax0)
	cwtmatran = signal.cwt(ranfirstDer, signal.ricker, widths)
	ax6.imshow(cwtmatran,cmap='RdPu',extent=[0, (GlobalVariables.num-GlobalVariables.window), 1, endRange],aspect='auto',vmax=abs(cwtmatran).max(), vmin=-abs(cwtmatran).max())
	ax6.set_xlabel('Position',size=6)
	ax6.set_yticks(ax5.get_yticks()[::2])
	ax6.set_title('Continuous Wavelet Transformation Convolved Over Range {0}-{1} for the First Derivative, Random Regions'.format(widths[0],endRange),size=8)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	print 'Plotted continuous wavelet transformation for elements and random regions'
	
	gs = gridspec.GridSpec(3,3,height_ratios=[2,1,1])
	gs.update(hspace=.65)
	
	# Short Fourier Transform
	ax7 = plt.subplot(gs[0,:],sharex=ax0)
	sbins = 30
	f1, t1, Zxx1 = signal.stft(firstDer,fs=1.0, window='hann',nperseg=sbins,noverlap=None)#,nperseg=11,noverlap=5
	ax7.pcolormesh(t1,f1,np.abs(Zxx1),cmap='RdPu')
	ax7.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#5fc85b')
	ax7.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#5fc85b')
	ax7.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#96c85b')
	ax7.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#96c85b')
	ax7.set_ylabel('Frequency',size=8)
	ax7.set_xlabel('Position',size=6)
	ax7.set_yticks(ax7.get_yticks()[::2])
	ax7.set_title('Short Fourier Transform over {0} bins for Elements'.format(sbins),size=8)
	
	# First Derivative
	ax8 = plt.subplot(gs[1,:],sharex=ax0)
	ax8.plot(GlobalVariables.fillX,firstDer,linewidth=1)
	ax8.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax8.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax8.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax8.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax8.axvspan(GlobalVariables.window,GlobalVariables.plotLineLocationThree,label='',alpha=0.1,facecolor = '#863eae')
	ax8.axvspan(GlobalVariables.plotLineLocationThree,GlobalVariables.plotLineLocationFour,label='',alpha=0.1,facecolor = '#ae3e9e')
	ax8.axvspan(GlobalVariables.plotLineLocationFour,(GlobalVariables.num-GlobalVariables.window-GlobalVariables.window),label='',alpha=0.1,facecolor = '#ae3e66')
	ax8.set_yticks(ax8.get_yticks()[::2])
	ax8.set_xlabel('Position',size=6)
	ax8.set_ylabel('Amplitude',size=8)
	ax8.set_title('First Derivative of Fitted Mean for Elements',size=8)
	
	ysdElement = get_element_coordinates(firstDer)
	frq2sd, Y2sd = perform_fourier_transfroms(ysdElement)
	ysdUp = get_up_stream_flank_coordinates(firstDer)
	frq3sd, Y3sd = perform_fourier_transfroms(ysdUp)
	ysdDown = get_down_stream_flank_coordinates(firstDer)
	frq4sd, Y4sd = perform_fourier_transfroms(ysdDown)
	print 'Performed fourier transform for elements'
	
	#FFT for sections of the smoothed second derivative
	ax9 = plt.subplot(gs[2,0])
	ax9.plot(frq3sd,abs(Y3sd),linewidth=1, color='#863eae')
	ax9.set_ylabel('|Y(freq)|',size=8)
	ax9.set_xlabel('Freq(Hz)',size=6) #AT Rate Change
	ax9.set_yticks(ax9.get_yticks()[::2])
	ax10 = plt.subplot(gs[2,1],sharey=ax9)
	plt.setp(ax10.get_yticklabels(), visible=False)
	ax10.plot(frq2sd,abs(Y2sd),linewidth=1, color='#ae3e9e')
	ax10.set_title('Power Series for Highlighted Regions, Elements',size=8)# Power Spectrum Analysis for FFT
	ax10.set_xlabel('Freq(Hz)',size=6)
	ax11 = plt.subplot(gs[2,2],sharey=ax9)
	plt.setp(ax11.get_yticklabels(), visible=False)
	ax11.plot(frq4sd,abs(Y4sd),linewidth=1, color='#ae3e66')
	ax11.set_xlabel('Freq(Hz)',size=6)
	
	sns.despine()
	plt.savefig(pp, format='pdf')
	print 'Plotted short fourier transform and fast fourier transform for elements'
	
	gs = gridspec.GridSpec(3,3,height_ratios=[2,1,1])
	gs.update(hspace=.65)

	ax12 = plt.subplot(gs[0,:],sharex=ax0)
	sbins = 30
	ranf1, rant1, ranZxx1 = signal.stft(ranfirstDer,fs=1.0, window='hann',nperseg=sbins,noverlap=None)#,nperseg=11,noverlap=5
	ax12.pcolormesh(rant1,ranf1,np.abs(ranZxx1),cmap='RdPu')
	ax12.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#5fc85b')
	ax12.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#5fc85b')
	ax12.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#96c85b')
	ax12.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#96c85b')
	ax12.set_ylabel('Frequency',size=8)
	ax12.set_xlabel('Position',size=6)
	ax12.set_yticks(ax12.get_yticks()[::2])
	ax12.set_title('Short Fourier Transform over {0} bins for Random Regions'.format(sbins),size=8)
	
	# First Derivative
	ax13 = plt.subplot(gs[1,:],sharex=ax0)
	ax13.plot(GlobalVariables.fillX,ranfirstDer,linewidth=1)
	ax13.axvline(x=GlobalVariables.plotLineLocationOne,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax13.axvline(x=GlobalVariables.plotLineLocationTwo,linewidth=.05,linestyle='dashed',color='#e7298a')
	ax13.axvline(x=GlobalVariables.plotLineLocationThree,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax13.axvline(x=GlobalVariables.plotLineLocationFour,linewidth=.05,linestyle='dashed',color='#bd4973')
	ax13.axvspan(GlobalVariables.window,GlobalVariables.plotLineLocationThree,label='',alpha=0.1,facecolor = '#863eae')
	ax13.axvspan(GlobalVariables.plotLineLocationThree,GlobalVariables.plotLineLocationFour,label='',alpha=0.1,facecolor = '#ae3e9e')
	ax13.axvspan(GlobalVariables.plotLineLocationFour,(GlobalVariables.num-GlobalVariables.window-GlobalVariables.window),label='',alpha=0.1,facecolor = '#ae3e66')
	ax13.set_yticks(ax13.get_yticks()[::2])
	ax13.set_xlabel('Position',size=6)
	ax13.set_ylabel('Amplitude',size=8)
	ax13.set_title('First Derivative of Fitted Mean for Random Regions',size=8)
	
	ranysdElement = get_element_coordinates(ranfirstDer)
	ranfrq2sd,ranY2sd = perform_fourier_transfroms(ranysdElement)
	ranysdUp = get_up_stream_flank_coordinates(ranfirstDer)
	ranfrq3sd,ranY3sd = perform_fourier_transfroms(ranysdUp)
	ranysdDown = get_down_stream_flank_coordinates(ranfirstDer)
	ranfrq4sd,ranY4sd = perform_fourier_transfroms(ranysdDown)
	print 'Performed fourier transform for random regions'

	#FFT for sections of the smoothed second derivative
	ax14 = plt.subplot(gs[2,0])
	ax14.plot(ranfrq3sd,abs(ranY3sd),linewidth=1, color='#863eae')
	ax14.set_ylabel('|Y(freq)|',size=8)
	ax14.set_xlabel('Freq(Hz)',size=6) #AT Rate Change
	ax14.set_yticks(ax14.get_yticks()[::2])
	ax15 = plt.subplot(gs[2,1],sharey=ax9)
	plt.setp(ax15.get_yticklabels(), visible=False)
	ax15.plot(ranfrq2sd,abs(ranY2sd),linewidth=1, color='#ae3e9e')
	ax15.set_title('Power Series for Highlighted Regions, Random Regions',size=8)
	ax15.set_xlabel('Freq(Hz)',size=6)
	ax16 = plt.subplot(gs[2,2],sharey=ax9)
	plt.setp(ax16.get_yticklabels(), visible=False)
	ax16.plot(ranfrq4sd,abs(ranY4sd),linewidth=1, color='#ae3e66')
	ax16.set_xlabel('Freq(Hz)',size=6)

	sns.despine()
	pp.savefig()
	pp.close()
	print 'Plotted short fourier transform and fast fourier transform for random regions'

def main(dfWindow,names,ranWindow,fileName):
	print 'Running graph_signal_linesLibrary'
	ATgroup = graph_signal_lines(dfWindow,names,ranWindow,fileName)

if __name__ == "__main__":
	main()
