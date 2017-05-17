#Type something like this from within python: ARGS="2015_04_07_02 M 30 none 10.5 DM"
#template: ARGS="filename sex age drugs theta region"
#then type: execfile('PopSpikeAnal.py')
#If not in python, type: python PopSpikeAnal.py 031011_5Or M 30 none 10.5 DM
import os

#if running from windows, might need to execute the following command:
#os.getcwd()
#might need to use forward slashes in this command:
#pythondir="C:\Users\avrama\Documents\StochDiffcode\Python"
#os.chdir(pythondir)

import numpy as np
from matplotlib import pyplot
import sys
from pprint import pprint as pp
import glob
import argparse
import pickle
from scipy import optimize
import pop_spike_utilities as psu

user="VL"
#TWEAK THESE as NEEDED: If mislabeling Features
#1. such as labeling FV as Pop spike, make FVwidth larger
#2. labeling the end of the artifact as the popspike (make artdecaytime larger)
# Units are msec, mV
artdecaytime=2.0 #units: ms
FVwidth=0.5 #units: ms
noisethresh=1.5 #units: mV, if pospeak is this much greater than baseline, there is a problem
artifactthreshold=1.5 #units: mV make smaller if artifact is too small
big_popspike_factor=3  #possible problem if popspike is bigger than mean + factor*stdev
slope_std_factor=2 #<<<<<<<<<<<<<<<<<<< if slope +/- slope_std_factor*std excludes 0,  print warning

#####parameters that you may want to tweak ###############
#If induction is not during the last pause in recording, you need to fix some of the code below
first_peak_end_fraction=0.725      #do not look past this fraction of trace for the 1st popspike (maybe there is a second one that we will look for later).
artifact_window=0.5     #if artifact is not found in first artifact_window of trace, there is problem
baseline_minutes=15
tracesPerMinute=2      #how often you stimulate and sample per minute

# for some summary measures, average over a x min window (where x = sample window)
#surronding the specified sample_time postinduction
#Analyze these measures using ANOVA
sample_window=5
sample_times=[30,60,90,120]

if user=="AK":
    datadir="F:\AlexData//"
    #where to put pickel files.  Make sure this directory exists or the program won't work.
    outputdir="F:\AlexData/pickle/"
if user=="VL":
    #Directory where data is located
    datadir="E:\FIELDS/ValerieData//"
    #where to put pickel files, relative to current path.  Make sure this directory exists or the program won't work.
    outputdir="C:\Users/Sarah/My Documents/Python Scripts/Pickle//"

#plasticity induction produces a time gap in the time wave, more than 10% greater than time between traces
induction_gap=(60./tracesPerMinute)*1.1       #10% longer than the seconds between popspikes
pretraces=tracesPerMinute*baseline_minutes
two_hours=120
sec_per_msec=0.001
baselinepoints=10	#This is just for finding the artifact
plotYN=1

###################### below here, the only parameter is experiment naming convention #######################
anal_params={'artdecay':artdecaytime, 'FVwidth':FVwidth, 'noisethresh':noisethresh, 'big_popspike_factor': big_popspike_factor,
             'artifactthresh': artifactthreshold,'artifact_wind':artifact_window,'1st_peak_end': first_peak_end_fraction,
             'induction_gap': induction_gap, 'baseline_min': baseline_minutes, "baselinepts": baselinepoints}


try:
	commandline = ARGS.split() #in python: define space-separated ARGS string
	do_exit = False
except NameError: #if you are not in python, read in filename and other parameters
	commandline = sys.argv[1:]
	do_exit = True

params=psu.parse_args(commandline,do_exit,0)
	
#Looks for file specified experiment name
filenamepattern=datadir+params.exper+".lvm"
print "Looking for files using: ", filenamepattern
filename = glob.glob(filenamepattern)
if (len(filename)==0):
	print "You mistyped the filenames. Python found no such files:"
	pp(filename)

#read in data file
with open(filename[0],'r') as f:
    for line in f:
        if "Samples" in line:
            break
timepoints_per_trace=int(line.split()[1])
first_peak_end=int(first_peak_end_fraction*timepoints_per_trace)
print "timepoints_per_trace",timepoints_per_trace
latest_artifact=artifact_window*timepoints_per_trace 
data = np.loadtxt(filename[0], skiprows=24) 

#first column of datafile is time, second column is Vm
wholetime=data[:,0]
tempVm=data[:,1]

dt=wholetime[1]-wholetime[0]
artifactdecay=int((artdecaytime*sec_per_msec)/dt)
if FVwidth>0:
        FVwidth_pts=int((FVwidth*sec_per_msec)/dt)
else:
        FVwidth_pts=0

#reshape the data to put each trace in separate column
datalength=np.shape(tempVm)[0]
numtraces=datalength/timepoints_per_trace
Vm_traces=np.reshape(tempVm, (numtraces, -1))



#create time wave that increments by dt, for only 400 samples
time=wholetime[0:timepoints_per_trace]

#extract the start time for each trace
tracetime=np.zeros(numtraces)
for i,j in enumerate(range(0, datalength, timepoints_per_trace)):
        tracetime[i]=wholetime[j]

#find the induction time - which is the trace after the first large gap in recordings
#what if there is more than one gap?  The last gap is induction gap
trace_diff=np.diff(tracetime)
gaps=np.where(trace_diff>induction_gap)[0]
print "gaps in recordings", gaps
if len(gaps)==0:
        print "whoops, no gap, possibly no induction"
else:
        induction=gaps[-1]

#if gap is position/array index n, then the gap is between n and n+1 in the tracetime array
#column n is last trace of baseline, column n+1 contains first trace of follow-on
#To analyze baseline, uses traces from n-(baseline_minutes*2) to n 
#To analyze after induction, use traces from n+1 to the end

if induction<pretraces:
        pretraces=induction
starttrace=induction-pretraces+1
num_usable_traces=min(numtraces-starttrace,(two_hours+baseline_minutes)*tracesPerMinute)
num_usable_minutes=int(np.ceil(float(num_usable_traces)/tracesPerMinute))
goodtraces=range(starttrace,starttrace+num_usable_traces)
print "traces",numtraces,"usable", num_usable_traces,"start", starttrace,"two hours", (two_hours+baseline_minutes)*tracesPerMinute,"last usable", goodtraces[-1]

#initialize some arrays to store trace characteristics
base=np.zeros((num_usable_traces))
peak=np.zeros((num_usable_traces))
peaktime=np.zeros((num_usable_traces))
pospeak=np.zeros((num_usable_traces))
pospeaktime=np.zeros((num_usable_traces))
FVsize=np.zeros((num_usable_traces))
amp=np.zeros((num_usable_traces))
baseline_start=np.zeros((num_usable_traces))
popspikestart=np.zeros((num_usable_traces))

#Analyze all traces from start to end, start is induction traces prior to gap
#pyplot.ion()
fig,axes=pyplot.subplots()
fig.canvas.set_window_title('Problem Traces for '+params.exper)
axes.set_ylabel('Vm (mV)')
axes.set_xlabel('Time (sec), red star marks end of artifact decay')
problem=0
for tracenum in goodtraces:
    index=tracenum-starttrace
    #Find artifact, where trace exceeds tempbaseline by some threshold
    artifactbegin=0
    tempbaseline=np.mean(Vm_traces[tracenum,0:baselinepoints])
    over_thresh=np.abs(Vm_traces[tracenum]-tempbaseline)>artifactthreshold
    #over_thresh=np.where(Vm_traces[tracenum]-tempbaseline>artifactthreshold)
    if(np.all(-over_thresh)):
	prob=1
    else:
	prob=0
	artifactbegin=np.min(np.where(over_thresh))
	baseline_start[index]=artifactbegin-int(sec_per_msec/dt)
    	ps_start=int(artifactbegin+artifactdecay)
    	popspikestart[index]=ps_start*dt
    	#print(artifactbegin, latest_artifact, baseline_start[index], ps_start, popspikestart[index])
    if artifactbegin==0 or artifactbegin>latest_artifact or prob:
        print "WARNING!!!!! NO ARTIFACT FOUND for trace",tracenum,'artifact', artifactbegin*dt, Vm_traces[tracenum,artifactbegin]-tempbaseline
        problem+=1
        base[index]=np.nan
        peaktime[index]=np.nan
        peak[index]=np.nan
        pospeaktime[index]=np.nan
        pospeak[index]=np.nan
        amp[index]=np.nan
        FVsize[index]=np.nan
        if problem<10:
            axes.plot(time,Vm_traces[tracenum,:]+index*0.1,label=tracenum)
            axes.plot(popspikestart[index],Vm_traces[tracenum,ps_start]+index*0.1,'r*')
            fig.canvas.draw()
    else:
        #Artifact found, start analyzing data 1 ms prior to artifact
        #Calculate base as mean of 1 ms prior to artifact
        base[index]=np.mean(Vm_traces[tracenum,baseline_start[index]:artifactbegin])
        #find negative peak and peak location, beginning at artifactdecay+FVwidth
        peakpoint=Vm_traces[tracenum,ps_start+FVwidth_pts:first_peak_end].argmin()+ps_start+FVwidth_pts
        peaktime[index]=peakpoint*dt
        peak[index]=Vm_traces[tracenum,peakpoint]
        #find the peak which divides fiber volley and popspike
        #if the negative peak is found as part of artifact decay, plot as problem trace
        if (peakpoint-ps_start)<2:
            print "WARNING!!!!! peakpoint found at end of artifact decay for trace",tracenum,'artifact:', artifactbegin*dt, 'peak:', peakpoint*dt 
            axes.plot(time,Vm_traces[tracenum,:]+index*0.1,label=tracenum)
            axes.plot(popspikestart[index],Vm_traces[tracenum,ps_start]+index*0.1,'r*')
            axes.plot(popspikestart[index],Vm_traces[tracenum,peakpoint]+index*0.1,'r*')
            fig.canvas.draw()
            pospeaktime[index]=np.nan
            pospeak[index]=np.nan
            amp[index]=np.nan
            FVsize[index]=np.nan
            problem+=1
        else:
            pospeakpoint=Vm_traces[tracenum,ps_start:peakpoint].argmax()+ps_start
            #print tracenum,pospeakpoint
            pospeaktime[index]=pospeakpoint*dt
            pospeak[index]=Vm_traces[tracenum,pospeakpoint]
            #find fiber volley size - predefined width
            FVsize[index]=np.abs(Vm_traces[tracenum,ps_start:ps_start+FVwidth_pts].min()-base[index])
	    #If pospeak is large noise artifact, this is problem trace
            if pospeak[index]-base[index]>noisethresh:
                    print "WARNING!!!! positive peak is really big.  Is this OK?"
                    axes.plot(time,Vm_traces[tracenum,:]+index*0.1,label=tracenum)
                    axes.plot(pospeaktime[index],Vm_traces[tracenum,pospeakpoint]+index*0.1,'mo')
                    axes.plot(peaktime[index],Vm_traces[tracenum,peakpoint]+index*0.1,'g*')
                    pospeaktime[index]=np.nan
                    pospeak[index]=np.nan
                    amp[index]=np.nan
                    FVsize[index]=np.nan
                    problem+=1
            else:
                    #calculate amplitude as difference between negative peak and either baseline or positive peak
                    if base[index] > pospeak[index]:
                            amp[index]=(base[index]-peak[index])
                    else: 
                            amp[index] = (pospeak[index]-peak[index])
axes.legend(fontsize=8, loc='best')

meanpre=amp[~np.isnan(amp[0:pretraces])].mean()
pretraces=baseline_minutes*tracesPerMinute
for index in range(pretraces,len(goodtraces)):
    #find nanmax (and location of nanmax)
    #calculate mean and std without the max, compare max to criteria, but using stdev =4
    #if too high, exclude, repeat until the answer is no.
    criteria=np.nanmean(amp[pretraces:]/meanpre)+big_popspike_factor*np.nanstd(amp[pretraces:]/meanpre)
    if np.abs(amp[index]/meanpre)>np.abs(criteria):
	problem+=1
	print "Possible PROBLEM with amplitude of peak for trace", goodtraces[index], "index", index, "normamp",amp[index]/meanpre,"crit",criteria
	amp[index]=np.nan
        st=baseline_start[index]
        end=int((peaktime[index]+0.004)/dt)
        axes.plot(time[st:end],Vm_traces[goodtraces[index],st:end],'-',label=tracenum)
#if problem==0:
#    pyplot.close(fig)
#else:
#    text=raw_input('ready for good traces plot? (y/n)')
#    pyplot.close(fig)

#plot traces with extracted peaks to verify
if plotYN:
        fig,axes=psu.plot_peaks(params.exper,time,Vm_traces,peaktime,pospeaktime,peak,pospeak,popspikestart,base,goodtraces,baseline_start,FVwidth*sec_per_msec)

#Next, normalize amplitude by mean popspike from preinduction traces, and
# average two samples per minute to produce a single value per minute, both for popspike and fibre volley
popspikenorm=np.zeros((num_usable_minutes))
FVnorm=np.zeros((num_usable_minutes))
popspikeminutes=np.zeros((num_usable_minutes))
meanFVpre=FVsize[0:pretraces].mean()
for sample in range(num_usable_minutes):
	popspikenorm[sample]=np.nanmean(amp[sample*tracesPerMinute:(sample+1)*tracesPerMinute])/meanpre	
	popspikeminutes[sample]=(tracetime[starttrace+sample*tracesPerMinute]-tracetime[induction])/60        #60 sec per minute
	if FVwidth > 0:
		FVnorm[sample]=FVsize[sample*tracesPerMinute:(sample+1)*tracesPerMinute].mean()/meanFVpre
	else:
		FVnorm[sample]=0

#Extract a few measures for statistical analysis, e.g. extract the 30, 60, 90, 120 min PopSpike means and FiberVolley
#Times are specified in sample_times array 
popspike_timesamples=np.zeros((len(sample_times)+1))
FV_timesamples=np.zeros((len(sample_times)+1))
popspike_timesamples[0]=meanpre
FV_timesamples[0]=meanFVpre
deletetime=[]
for index in range(0,len(sample_times)):
        norm_index=sample_times[index]+baseline_minutes
        if len(popspikenorm)<norm_index:
                print "calculating time samples; recording of ", popspikeminutes[-1], "min shorter than ", sample_times[index]
                popspike_timesamples[index+1]=np.nan
                FV_timesamples[index+1]=np.nan
                deletetime.append(index)
        else:
                popspike_timesamples[index+1]=np.nanmean(popspikenorm[norm_index-sample_window:norm_index])
                FV_timesamples[index+1]=FVnorm[norm_index-sample_window:norm_index].mean()
#May need to delete this line if column_stack in GrpAvgPopSpike complains about different length arrays
new_sample_times=np.delete(sample_times,deletetime)+(baseline_minutes-sample_window/2)

print "ready for summary plot, and fit line to baseline"
#Plot popspikenorm vs minutes, fit line to baseline, print slope
validbasepts=~np.isnan(popspikenorm[0:baseline_minutes])
popt,pcov=optimize.curve_fit(psu.line,popspikeminutes[validbasepts],popspikenorm[validbasepts])
Aopt,Bopt=popt
Astd,Bstd=np.sqrt(np.diag(pcov))
if  np.abs(Bopt)-slope_std_factor*Bstd>0:
        print "********************** WARNING ************** \n Baseline Slope of ", Bopt,"+/- 2*",Bstd, "sig diff than zero"

if plotYN:
        psu.plot_summary(popspikeminutes,popspikenorm,baseline_minutes,popspike_timesamples,new_sample_times,Aopt,Bopt,fig,axes)

#tracedict includes all traces that you potentially want to analyze by groups in GroupAvg
#modify argsdict to include all the other bits of information you want, and change dictionary keys to be meaningful for you.
tracedict={'amp':amp,'peaktime':peaktime,'FVsize':FVsize,'popspikenorm':popspikenorm,'popspikeminutes':popspikeminutes, 'FVnorm': FVnorm,'PS_mean':popspike_timesamples,'FV_means':FV_timesamples,'slope':Bopt,'slope_std':Bstd}

#This writes the file for GrpAvg.py
#it includes all the experiment parameters
datadict = dict(trace=tracedict, #datadict name wont be accessible outside this program but keys are
		parameters=params, anal_params=anal_params) #first is key second is value
outpickle = outputdir+params.exper + '.pickle'
with open(outpickle, 'w') as f:
	pickle.dump(datadict, f) #saves the contents of datadict in file "f"... 

if plotYN:
        # show the graphs
        pyplot.show()
