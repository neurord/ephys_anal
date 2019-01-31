'''
Correlogram.py: calculates shuffle corrected cross-correlogram, averaged across all neuron pairs
Usage: python3 Correlogram.py --file filename
Optional parameters specifying start and end time, e.g.
python3 Correlogram.py --file Ethanol10min.csv --tstart 19.792 --tend 79.792
if start and end time of spikes to be analyzed are the same for several files, can process multiple files at once:
python3 Correlogram.py --file Ethanol10min.csv --file Saline10min.cvs
Additional optional parameters: number of histogram bins, and maximum bin delta T. Default is binMax=50e-3
fileformat:
if csv file, each column is the spike time of a different neuron. Software detects how many neurons
alternatively, all data in two columns, with 1st column the spike time, and 2nd column the number.  I
default number of shuffles=5
This does not return correlation between each pair of neurons
   between a neuron and all others, implement by not summing over the neurons
   among low firing or high firing: implement by creating two subsets and processing separately
   between pairs of neurons: probably need to eliminate broadcast and calculate time differences using loops
      e.g., calculate dt and histogram between pairs of neurons.  Then sum over all pairs or subset of pairs
   isi_histogram_by_number_plot seems to calculate histogram for individual spikes, BUT, may not be correct and needs weight - this is good starting point.
'''
from __future__ import print_function as _, division as _
import sys
import os
import numpy as np

#from concurrent import futures
import argparse
import shuffle
import histograms as hist
import plot_utils as pu
import util

timetype = 'f4'
histtype = 'u4'

plotspikes=0
decimals=6
#default parameters for doing the histogram
nBinsHist=100
binMax=50e-3
binMin=-binMax

#parameters for unused shuffle methods
isiBurstCutoff = 0.025
firstSpikeShift = 0.015
reshuffleOffset = 0.5

#parameters
#number of times to shuffle data for shuffle subtraction
shuffles = 5
#ignore spikes that occur prior to network equilibration
steadyState = 0.0 # 0.3
#only analyze spikes that occured prior to the network equilibration
spikedBeforeSteadyOnly = False
#only use neurons with at least this many spikes
minSpikeCount = 0 # 3

#parameters used for histogram methods
useFirstInHistogram = True
useLastInHistogram = True

#number of neurons for test data
#input parameters
totalTime=1000 #maximum time stamp of spikes
netsize=1000

#Note: saline: 11.244 to 71.244
#ethanol from 19.792 to 79.792

parser = argparse.ArgumentParser()
parser.add_argument('--file', nargs='*') 
parser.add_argument('--tstart', default=0,type=float,help='starting time of usable spike times')
parser.add_argument('--tend', default=totalTime,help='end time of usable spike times',type=float)
parser.add_argument('--binMax', default=binMax,help='maximum delta T for histogram',type=float)
parser.add_argument('--nBins', default=nBinsHist,help='number of bins for histogram',type=int)

try:
    args=ARGS.split(" ")
    print ("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

p=parser.parse_args(args)
print('files to process:', p.file)
binWidth=(binMax-binMin)/nBinsHist
#print "binWidth:",binWidth
bins=np.linspace(binMin, binMax, nBinsHist)
filename=p.file[0]

if filename == 'flat':
    print('flat')
    maxspikes = 5
    spikeTime = [sorted(np.random.rand(np.random.randint(maxspikes)) * binMax)
                 for _ in range(netsize)]
    filename = 'flat'

elif filename == 'flat_with_cutoff':
    meanspikes = 1200
    min_isi=0.004
    spikeTime = [sorted(np.random.rand(np.random.poisson(meanspikes)) * binMax)
                 for _ in range(netsize)]
    for i, tt in enumerate(spikeTime):
        tt = np.array(spikeTime[i])
        spikeTime[i] = [tt[0]] + list(tt[1:][np.diff(tt) > 0.004])

elif filename == 'flat_with_trains':
    spikeTime = [sorted(np.hstack(util.trains()))
                 for _ in range(netsize)]

else:
    filetype=filename.split('.')[-1]
    if filetype=='csv': #one column per neuron
        spikeTime=[]
        netsize=0
        for fname in p.file:
            netdata = util.readcsv(fname)
            subnetsize=np.shape(netdata)[1]
            netsize=netsize+subnetsize
            #next, create one list per neuron, exclude spike times greater than 1 minute
            for neurnum in range(subnetsize):
                neur_spikes=[st for st in netdata[:,neurnum] if (st<p.tend and st > p.tstart)]
                spikeTime.append(neur_spikes)
        print (fname,'netsize',np.shape(netdata)[1], 'total netsize', netsize, np.shape(spikeTime)) 
    else: #2 columns of data: col 1: time of spike, col 0: neuron number of spike
        datas=[]
        netsizes=[0]
        for fname in p.file:
            data=np.loadtxt(fname)
            datas.append(data)
            netsizes.append(int(np.max(data[:,0])-np.min(data[:,0]))+1)
        netsize=np.sum(netsizes)
        spikeTime=[[] for _ in range(netsize)]
        for ns,neurdata in enumerate(datas):
            for i in range(neurdata.shape[0]):
                spikeTime[int(neurdata[i,0]) - 1+netsizes[ns]].append(neurdata[i,1])
    # #Removing initial delay from entire dataset to prevent initial offset from biasing shuffling
    # baseline = data[0,1]
    # for i in range(netsize):
    #     for j in range(len(spikeTime[i])):
    #         spikeTime[i][j] = round(spikeTime[i][j]-baseline,decimals)

    #del data

########### create array from lists, use np.nan to fill array for neurons with fewer spikes
T = np.empty((netsize, max(len(l) for l in spikeTime)), dtype=timetype)
for i, l in enumerate(spikeTime):
    T[i, :len(l)] = l
    T[i, len(l):] = np.nan
#if equilibration time is defined, subtract it from all spikeTimes; don't use initial spikeTime if negative
if steadyState:
    T -= steadyState
    #This next replaces everything with nans, except spikes earlier than equilibrium time
    if spikedBeforeSteadyOnly:
        T[np.nanmin(T, axis=1) > 0] = np.nan
    #don't use negative spikeTimes
    T[T <= 0] = np.nan
    T = np.sort(T, axis=1)

if minSpikeCount > -1:
    #exclude neurons whose number of spikes is less than the minimum
    T = T[ ~np.isnan(T[:, minSpikeCount])]
    netsize = T.shape[0]

# Shuffle into newT
print('begin shuffle')
newT = shuffle.shuffle_everything(T,shuffles,netsize)
#newT = shuffle.shuffle_bursts(T,shuffles,netsize,isiBurstCutoff)
#newT = shuffle.shuffle_everything_first_long(T,shuffles,netsize,isiBurstCutoff)
#newT = shuffle.shuffle_everything_keep_initial(T,shuffles,netsize)
#newT = shuffle.shuffle_everything_keep_initial_gaussian_shift(T,shuffles,netsize, firstSpikeShift)
#newT = shuffle.shuffle_everything_twice(T,shuffles,netsize, reshuffleOffset)

numspikes = (~np.isnan(T)).sum(axis=1)
print('numspikes', numspikes)
#weight is 1/(N_A*N_B), doesn't include 2*deltaT
weights = 1 / (numspikes[:,None] * numspikes)
weights[np.isinf(weights)] = 0
spikers = shuffle.nonnancount(T)

#hist determines number of coincident spikes within a X ms - X is the bin
hists=[hist.diff_and_hist3((newT, i, bins, weights, useFirstInHistogram, useLastInHistogram) ) for i in range(netsize)]
hist_shuffle = np.sum((h[None, :] for h in hists), dtype=float, axis=1)[0] / netsize

#doing calculations in parallel doesn't work well with large number of spikes:
'''
number_cpus = os.sysconf('SC_NPROCESSORS_ONLN')
with futures.ThreadPoolExecutor(number_cpus) as exe:
    hists = exe.map(hist.diff_and_hist3,
                    ((newT, i, bins, weights, useFirstInHistogram, useLastInHistogram) for i in range(netsize)))
    hist_shuffle = np.sum((h[None, :] for h in hists), dtype=float, axis=1)[0] / netsize
    del hists
'''
hist_s = hist_shuffle/shuffles
print('shuffle',hist_s)

oldT = T[None, ...]
hists=[hist.diff_and_hist3((oldT, i, bins, weights, useFirstInHistogram, useLastInHistogram) ) for i in range(netsize)]
hist_unshuffled = np.sum((h[None, :] for h in hists), dtype=float, axis=1)[0] / netsize
'''
oldT = T[None, ...]
with futures.ThreadPoolExecutor(number_cpus) as exe:
    hists = exe.map(hist.diff_and_hist3,
                    ((oldT, i, bins, weights, useFirstInHistogram, useLastInHistogram) for i in range(netsize)))
    hist_unshuffled = np.sum((h[None, :] for h in hists), dtype=float, axis=1)[0] / netsize
    del hists
'''
print('data',hist_unshuffled)

hist_diff = hist_unshuffled - hist_s
middles = (bins[:-1] + bins[1:]) / 2
#better determine percentile from all 5 individual shuffled histograms
tile95=np.percentile(np.sort(hist_s),95)
sig=np.where(hist_diff>tile95)
sig_bin=middles[sig]
sig_val=hist_diff[sig]
print("SIGNIFICANT VALUES, 95%", tile95)
for b,v in zip(sig_bin,sig_val):
    print(b,':',v)

if True:
    pu.isi_histogram_plot(middles,  filename, shuffled=hist_s, unshuffled=hist_unshuffled)
    pu.isi_histogram_diff_plot(middles,  filename, hist_diff,tile95)
    pu.spike_ini_plot(T, filename, spikers)
    #pu.isi_histogram_by_number_plot(bins, T,  filename, additive=False)
    #pu.isi_recurrence_plot2(T, filename)
    #raw_input function seems to have disappeared
    #raw_input()

#save histograms as well as bins for plotting later
np.savetxt(os.path.basename(fname)+"_hist.txt",hist_diff)

mouse=os.path.basename(fname).split('.')[0].split('_')[0]
np.savetxt(mouse+'_bins.txt', middles)
#

