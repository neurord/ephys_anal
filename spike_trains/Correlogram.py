'''
Correlogram.py: calculates shuffle corrected cross-correlogram, averaged across all neuron pairs
Usage: python3 Correlogram.py --file filename
Optional parameters specifying start and end time, e.g.
python3 Correlogram.py --file Ethanol10min.csv --tstart 19.792 --tend 79.792
if start and end time of spikes to be analyzed are the same for several files, can process multiple files at once:
python3 Correlogram.py --file Ethanol10min.csv --file Saline10min.cvs
fileformat:
if csv file, each column is the spike time of a different neuron. Software detects how many neurons
alternatively, all data in two columns, with 1st column the spike time, and 2nd column the number.  In this case must specify the number of neurons of each type.  Assumes that there of two neuron types.
default number of shuffles=5
min and max delta t for final graph specified as binMax
'''
from __future__ import print_function as _, division as _
import sys
import os
import numpy as np
import csv
from concurrent import futures
import argparse
import shuffle
import histograms as hist
import plot_utils as pu

timetype = 'f4'
histtype = 'u4'

plotspikes=0
decimals=6
#parameters for doing the histogram
nBinsHist=100
binMax=50e-3
binMin=-binMax
binWidth=(binMax-binMin)/nBinsHist
#print "binWidth:",binWidth
bins=np.linspace(binMin, binMax, nBinsHist)

#parameters for unused shuffle methods
isiBurstCutoff = 0.025
firstSpikeShift = 0.015
reshuffleOffset = 0.5

#parameters
#number of times to shuffle data for shuffle subtraction
shuffles = 5  
steadyState = 0.0 # 0.3
spikedBeforeSteadyOnly = False
#only use neurons with at least this many spikes
minSpikeCount = 0 # 3

#parameters used for histogram methods
useFirstInHistogram = True
useLastInHistogram = True

#number of neurons for test data
#input parameters
NumPairs=0
totalTime=1000
netsize=1000

#Note: saline: 11.244 to 71.244
#ethanol from 19.792 to 79.792

parser = argparse.ArgumentParser()
parser.add_argument('--file', nargs='*') 
parser.add_argument('--tstart', default=0,type=float)
parser.add_argument('--tend', default=totalTime,type=float)

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

def readcsv(fname):
    data=[]
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            data.append(row)
    neurnames=[name for name in data[0] if len(name)>0]
    numrows=len(data[1:])
    numcols=len(neurnames)
    print('file,num_neurons',fname, numcols, 'names',neurnames)
    timestamp=np.empty((numrows,numcols))
    for index,row in enumerate(data[1:]):
        timestamp[index]=[float(row[j]) if len(row[j])>0 else np.nan for j in range(numcols) ]
    return timestamp

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
    def trains():
        meantrains = 3
        meanspikes = 1.5
        intertrainisi = 0.012
        intertrainisinoise = 0
        intertrainnoise = 0.001
        min_isi = 0.004

        time = 0
        while time < binMax:
            t = np.random.exponential(binMax / meantrains)
            time += t

            n = np.random.geometric(1 / meanspikes)
            itisi = intertrainisi + np.random.randn(n) * intertrainisinoise
            isis = np.cumsum(min_isi + np.random.exponential(itisi - min_isi, size=n))
            times = time + isis
            yield times[(times < binMax) & (times >= 0)]
            time = times[-1]

    spikeTime = [sorted(np.hstack(trains()))
                 for _ in range(netsize)]

else:
    filetype=filename.split('.')[-1]
    if filetype=='csv': #one column per neuron
        spikeTime=[]
        netsize=0
        for fname in p.file:
            netdata = readcsv(fname)
            subnetsize=np.shape(netdata)[1]
            netsize=netsize+subnetsize
            #next, create one list per neuron, exclude spike times greater than 1 minute
            for neurnum in range(subnetsize):
                neur_spikes=[st for st in netdata[:,neurnum] if (st<p.tend and st > p.tstart)]
                spikeTime.append(neur_spikes)
        print (fname,'netsize',np.shape(netdata)[1], 'total netsize', netsize, np.shape(spikeTime)) 
    else: #2 columns of data: col 1: time of spike, col 0: neuron number of spike, only 2 files allowed
        filename2 = args[1] if len(args) > 1 else filename.replace('_D1_', '_D2_')
        data = np.loadtxt(filename)
        data2 = np.loadtxt(filename2)
        netsize1=np.max(data[:,0])-np.min(data[:,0])
        netsize2=np.max(data2[:,0])-np.min(data2[:,0])
        netsize=netsize1+netsize2
        spikeTime=[[] for _ in range(netsize)]
        for i in range(data.shape[0]):
            spikeTime[int(data[i,0]) - 1].append(data[i,1])
        for i in range(data2.shape[0]):
            spikeTime[int(data2[i,0]) - 1 + netsize1].append(data2[i,1])

    # #Removing initial delay from entire dataset to prevent initial offset from biasing shuffling
    # baseline = data[0,1]
    # for i in range(netsize):
    #     for j in range(len(spikeTime[i])):
    #         spikeTime[i][j] = round(spikeTime[i][j]-baseline,decimals)

    #del data, data2

#Old algorithm:
# Calculating correlation
#for i in range(netsize):
#    freq.append(1/mean(np.diff(spikeTime[i])))
#    meanISI.append(mean(np.diff(spikeTime[i])))
#    stdISI.append(std(np.diff(spikeTime[i])))
#    for j in range(netsize):
#	if (i != j):
#		timeDiffs=list()
#		for m in spikeTime[i]:
#			for n in spikeTime[j]:
                #Only include time differences within bins, but it might not matter
#				if (abs(m-n)<binMax):
#					timeDiffs.append(m-n)

#create array from lists, use np.nan to fill array for neurons with fewer spikes
T = np.empty((netsize, max(len(l) for l in spikeTime)), dtype=timetype)
for i, l in enumerate(spikeTime):
    T[i, :len(l)] = l
    T[i, len(l):] = np.nan
if steadyState:
    T -= steadyState
    if spikedBeforeSteadyOnly:
        T[np.nanmin(T, axis=1) > 0] = np.nan

    T[T <= 0] = np.nan
    T = np.sort(T, axis=1)

if minSpikeCount > -1:
    #exclude neurons whose number of spikes is less than the minimum
    T = T[-np.isnan(T[:, minSpikeCount])]
    netsize = T.shape[0]

# Shuffle into newT
print('begin shuffle')
newT = shuffle.shuffle_everything(T,shuffles,netsize)
#newT = shuffle.shuffle_bursts(T,shuffles,netsize,isiBurstCutoff)
#newT = shuffle.shuffle_everything_first_long(T,shuffles,netsize,isiBurstCutoff)
#newT = shuffle.shuffle_everything_keep_initial(T,shuffles,netsize)
#newT = shuffle.shuffle_everything_keep_initial_gaussian_shift(T,shuffles,netsize, firstSpikeShift)
#newT = shuffle.shuffle_everything_twice(T,shuffles,netsize, reshuffleOffset)

numspikes = (-np.isnan(T)).sum(axis=1)
print('numspikes', numspikes)
weights = 1 / (numspikes[:,None] * numspikes)
weights[np.isinf(weights)] = 0
spikers = shuffle.nonnancount(T)

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
print('shuffle',hist_shuffle / shuffles)

hists=[hist.diff_and_hist3((newT, i, bins, weights, useFirstInHistogram, useLastInHistogram) ) for i in range(netsize)]
hist_unshuffled = np.sum((h[None, :] for h in hists), dtype=float, axis=1)[0] / netsize
'''oldT = T[None, ...]
with futures.ThreadPoolExecutor(number_cpus) as exe:
    hists = exe.map(hist.diff_and_hist3,
                    ((newT, i, bins, weights, useFirstInHistogram, useLastInHistogram) for i in range(netsize)))
    hist_unshuffled = np.sum((h[None, :] for h in hists), dtype=float, axis=1)[0] / netsize
    del hists
'''
print('data',hist_unshuffled)

if True:
    pu.isi_histogram_plot(bins,  filename, shuffled=hist_shuffle / shuffles, unshuffled=hist_unshuffled)
    pu.isi_histogram_diff_plot(bins,  filename, hist_unshuffled, hist_shuffle / shuffles)
    pu.spike_ini_plot(T, filename, spikers)
    pu.isi_recurrence_plot2(T, filename)
    #raw_input function seems to have disappeared
    #raw_input()

hist_s = hist_shuffle/shuffles
diff = hist_unshuffled - hist_s
np.savetxt(os.path.basename(fname)+"_hist.txt",diff)

