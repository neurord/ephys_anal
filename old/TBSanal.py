#This program reads through all analyzed experiments, extracts the dt between PSP and spike
#Run with execfile('TBSanal.py').  ARGS: name of file with list of experiments
#may need to adjust min_height to exclude artifacts or include smaller spikes
#make pause=0 to pause between display of detected spikes for experiments
import os
from igor import binarywave
import numpy as np
from matplotlib import pyplot
import sys
from pprint import pprint as pp
import glob
import pickle
from detect import detect_peaks

##########################################################################
#location of data waves (TBS) and pickle files (which have exper char)
rootdir="/home/avrama/EphysDataAnal/"
#rootdir="./"
wavedir=rootdir+"TBSdata/"
TBS_ending="_TBS_Waves/*"
stim_times=[0.0153,0.0353,0.0553,0.0753]
PSPdelay=0.004  #mean delay between E.C. stim and PSP onset
num_windows=10
pause=0
##########################################################################
#Function definitions
#Below puts your files from common exper in order by PGF and trace
def sortorder(fname):
    parts = fname.split('_')
    ans = int(parts[-3]), int(parts[-2])
    return ans

def setup(exper,numwindows):
    pyplot.ion()
    fig,axes=pyplot.subplots(numwindows,1)
    fig.canvas.set_window_title(exper)
    return fig,axes

def find_spikes(wave, min_height=-0.01, maximum_spike_rise_time=.004,threshval=0.02,min_rise=0.0003,refract=0.004):
    #wave[0]=time, wave[1] = Vm
    peaks = detect_peaks(wave[1])
    peaks = peaks[wave[1][peaks] > min_height]
    peaktime=wave[0][peaks]
    dt=wave[0][1]-wave[0][0]
    thresholds = np.empty(peaks.size)
    risetimes=np.empty(peaks.size)
    delete_list=[]
    for i in range(len(peaks)):
        start = max(0,(peaks[i] - maximum_spike_rise_time/dt))
        x = wave[0][start:peaks[i] + 1]
        y = wave[1][start:peaks[i] + 1]
        if len(y)>1:
            yderiv = np.diff(y)
            #spike threshold is point where derivative is 2% of steepest
            ythresh = threshval * yderiv.max()
            thresh= y[yderiv > ythresh].min()
            thresholds[i] = thresh
            threshtime=x[yderiv>ythresh].min()
            #calculate risetime for spike, and exclude if below min_rise (i.e., if artifact)
            risetime=peaktime[i]-threshtime
            if risetime<min_rise:
                #print "exclude: peaktime,risetime,peakval", peaktime[i],risetime,wave[1][peaks[i]]
                #exclude these stim artifact peaks from peaks
                delete_list.append(i)
        else:
           print "PROBLEM: i, numpeaks, peaktime",i, len(peaks),peaktime[i], "rise phase", start,"to", peaks[i]
           delete_list.append(i)
    if len(delete_list)>0:
        print "##artifact found, delete items:",delete_list,"from",peaks,
        peaks=np.delete(peaks,delete_list)
        thresholds=np.delete(thresholds,delete_list)
        peaktime=wave[0][peaks]
        print "    after",peaks
    #also need to exclude artifacts near end of spike being double counted, i.e. within refractory period
    spikefreq=np.diff(peaktime)
    #delete second peak, assuming this is the artifact.  Better method of determining which is artifact is needed
    dup=np.where(spikefreq<refract)[0]+1
    if len(dup)>0:
        print "-->dup found, delete items:",dup,"from",peaks,
        peaks=np.delete(peaks,dup)
        thresholds=np.delete(thresholds,delete_list)
        peaktime=wave[0][peaks]
        print "    after",peaks
    return peaks,thresholds

def exper_spikes(wavedir,expername,TBS_ending,axes,fig,stim_times):
    FileDir=wavedir+expername+TBS_ending
    TBSfiles=glob.glob(FileDir)
    print expername, "NUM FILES:", len(TBSfiles)
    if (len(TBSfiles)==0):
        print "No files in:", FileDir
        return 0,0
    else:
        TBSfiles = sorted(TBSfiles, key=sortorder)
        numspikes=np.zeros(len(TBSfiles))
        peaktime=[]
        #read in data file
        for i,filename in enumerate(TBSfiles):
            #read in data
            data = binarywave.load(filename)
            trace=data['wave']['wData']
            dt=data['wave']['wave_header']['hsA']
            stim_index=[int(round(st/dt)) for st in stim_times]
            npnts=data['wave']['wave_header']['npnts']
            endtime=len(trace)*dt
            tracetime=np.arange(0,endtime,dt) #array of points between 0 and endtime stepping by dt
            #
            #find peaks and peaktime using modification of Zbyszek's find_spikes
            peaks,thresholds=find_spikes([tracetime, trace])
            peaktime.append(tracetime[peaks])
            numspikes[i]=len(peaks)
            #plot trace and peaks to verify correct identification
            axes[(i/num_windows)].plot(tracetime,trace)
            for peak in peaks:
                if peak in stim_index:
                    axes[(i/num_windows)].plot(tracetime[peak],trace[peak],'^',label=filename.split('/')[-1].split('_')[5:7])
                    print filename.split('/')[-1], "peaks:",peaks,"check:", tracetime[peak]
                    axes[(i/num_windows)].legend(fontsize='xx-small')
                else:
                    axes[(i/num_windows)].plot(tracetime[peak],trace[peak],'*r')
        fig.canvas.draw()
        pyplot.show()
        return peaktime,numspikes

def deltat(spiketimes, PSPdelay, stim_times):
    min_interval=[]
    other_interval=[]
    all_intervals=[]
    for i in range(np.shape(spiketimes)[0]):
        min_interval.append([])
        other_interval.append([])
        all_intervals.append([])
        for j in range(np.shape(spiketimes[i])[0]):
            #positive delta means spiketime AFTER psp - forward pairing
            delta=spiketimes[i][j]-(np.array(stim_times)+PSPdelay)
            minloc=np.abs(delta).argmin()
            min_interval[i].append(delta[minloc])
            all_intervals[i].append(delta[minloc])
            #alternative value for spikes falling between two PSPs
            if spiketimes[i][j]>(stim_times[0]+PSPdelay) and spiketimes[i][j]<(stim_times[-1]+PSPdelay):
                if min_interval[i][j]>0:
                    other_interval[i].append(delta[minloc+1])
                    all_intervals[i].append(delta[minloc+1])
                else:
                    other_interval[i].append(delta[minloc-1])
                    all_intervals[i].append(delta[minloc-1])
            else:
                other_interval[i].append(min_interval[i][j])
    return min_interval,other_interval,all_intervals

def flat_stats(listoflists):
    flat_list=[x for y in listoflists for x in y]
    mean=np.mean(flat_list)
    stdev=np.std(flat_list)
    return mean,stdev

try:
	args = ARGS.split(",")
	print "ARGS =", ARGS, "commandline=", args
 	do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
	args = sys.argv[1:]
	print "commandline =", args
	do_exit = True

experlist=wavedir+args[0]      #E.g. PARAMSforSAS.txt
if str.find(args[0],"List")>0:
    outname=args[0][0:str.find(args[0],"List")]+'SAS.txt'
elif str.find(args[0],".txt")>0:
    outname=args[0][0:str.find(args[0],"List")]+'SAS.txt'
else:
    outname=args[0]+'SAS.txt'
print "outname=",outname
fig_spikes,axes_spikes=setup('allexpers',1)
fig,axes=setup('traces',num_windows)
dt_mean=[]
dt_std=[]
other_mean=[]
other_std=[]
all_mean=[]
all_std=[]
experout=[]
spikes_mean=[]
spikes_std=[]
#open experlist, read in line by line, skip the first line
for f in open(experlist,'r'):
    expername=f.split()[0]
    #
    [axes[i].clear() for i in range(len(axes))]
    fig.canvas.set_window_title(expername)
    #
    #read in data and extract spikes
    #
    spiketimes,numspikes=exper_spikes(wavedir,expername,TBS_ending, axes,fig,stim_times)
    experout.append(expername)
        #
    if sum(numspikes)!=0:
        #compare peaktimes with stimtime. 1: find minimum dt, 2: consider opposite for spikes b/t PSPs
        #
        min_interval,other_interval,all_intervals=deltat(spiketimes, PSPdelay, stim_times)
        #
        #Statistics (mean, stdev, histogram) of min_interval and other_interval
        #
        tmp_mean,tmp_std=flat_stats(min_interval)
        dt_mean.append(round(tmp_mean,8))
        dt_std.append(round(tmp_std,8))
        #tmp_mean,tmp_std=flat_stats(other_interval)
        #other_mean.append(round(tmp_mean,8))
        #other_std.append(round(tmp_std,8))
        tmp_mean,tmp_std=flat_stats(all_intervals)
        all_mean.append(round(tmp_mean,8))
        all_std.append(round(tmp_std,8))
    else:
        #if no spikes for that exper, append nans
        dt_mean.append(np.nan)
        dt_std.append(np.nan)
        all_mean.append(np.nan)
        all_std.append(np.nan)
    spikes_mean.append(numspikes.mean())
    spikes_std.append(numspikes.std())
    #
    #Graph at number of spikes over time, 1 trace per exper
    #
    axes_spikes.plot(range(len(numspikes)),numspikes)
    fig_spikes.canvas.draw()
    pyplot.show()
    if pause:
        text=raw_input('next? (y/n)')
        if text=='n':
            break
#
#ascii output of summary measures for statistical analysis
#
SASoutput=np.column_stack((experout,dt_mean,dt_std,spikes_mean,spikes_std,all_mean,all_std))
SASheader="exper              deltat_mean   deltat_std   spikes_mean   spikes_std   all_mean   all_std\n"
f=open(outname, 'w')
f.write(SASheader)
np.savetxt(f, SASoutput, fmt='%s', delimiter='   ')
f.close()
