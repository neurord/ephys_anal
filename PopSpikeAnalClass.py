#Type something like this from within python: ARGS="20160822orangedm105 M 30 none 10.5 DM"
#template: ARGS="filename sex age drugs theta region"
#then type:  
#If not in python, type: python PopSpikeAnal.py 031011_5Or M 30 none 10.5 DM
import os;
import numpy as np
from matplotlib import pyplot
import sys
from pprint import pprint as pp
import glob
import argparse
import pickle
from scipy import optimize
import pop_spike_utilities as psu
#unindent lines:

two_hours=120 #make global
sec_per_msec=0.001 #make global
sec_per_min=60.0 #make global  



class PopSpike:
   
    def __init__(self,params): #change to PopSpike(args)
        self.params = params
        self.sample_window = params.samp_win
        self.sample_times = params.samp_times
        self.filename_ending = params.file_end
        self.first_peak_end_fraction= params.first_peak_end_frac
        self.artifact_window =params.art_win
        self.baseline_minutes = params.base_min
        self.noisethresh= params.noisethres
        self.slope_std_factor=params.slope_std
        self.datadir = params.datadir
        self.outputdir = params.outputdir
        self.artdecaytime = params.decay*sec_per_msec
        self.artifactthreshold=params.artthresh
        self.FVwidth =params.FVwidth*sec_per_msec
        self.induction_gap_time=1.1 #multiplicative factor
        self.baselinepoints=10
        
        self.anal_params={'artdecay':self.artdecaytime, 'FVwidth':self.FVwidth, 'noisethresh':self.noisethresh, 
                 'artifactthresh': self.artifactthreshold,'artifact_wind':self.artifact_window,'1st_peak_end': self.first_peak_end_fraction,
                 'induction_gap': self.induction_gap_time, 'baseline_min': self.baseline_minutes, "baselinepts": self.baselinepoints}
        #10% longer than the seconds between popspikes       
        # arbitrary value, changed after running lookforfile  
        self.induction =0                
        self.problemtraces = {'tracenum' : [], 'artifactbegin' : [], 'tempbaseline' : []}
            
    def read_datafile(self):     
        filenamepattern=self.datadir+self.params.exper+self.filename_ending
        filename = glob.glob(filenamepattern)
        if (len(filename)==0):
            print ("You mistyped the filenames. Python found no such files:")
            pp(filename)           
            #read in data file. This is specific to labview output
        self.Vm_traces,self.tracetime,self.dt,self.time=psu.read_labview(filename[0])  
        self.tracesPerMinute=int(np.round(sec_per_min/(self.tracetime[-1]-self.tracetime[-2])))
        self.pretraces=self.tracesPerMinute*self.baseline_minutes
        self.induction_gap=(sec_per_min/self.tracesPerMinute)*1.1     
        self.numtraces=np.shape(self.Vm_traces)[0]
        self.timepoints_per_trace=np.shape(self.Vm_traces)[1]          
        print ("timepoints_per_trace",self.timepoints_per_trace, "Traces Per Minute", self.tracesPerMinute)                        
        text=input('OK to continue? (y/n)')
        if text=='n':
            print("exiting")
            raise Exception('wrong stimulation rate')
        #try to not use [index], do not have return from fuctions
        return
        
    def find_induction_gap(self):
        #convert_msec_to_points
        self.first_peak_end=int(self.first_peak_end_fraction*self.timepoints_per_trace)
        self.latest_artifact=self.artifact_window*self.timepoints_per_trace 
        self.artifactdecay=int(self.artdecaytime/self.dt)

        if self.FVwidth>0:
            self.FVwidth_pts=int(self.FVwidth/self.dt)
        else:
            self.FVwidth_pts=0                     

        self.trace_diff=np.diff(self.tracetime) #we dont use this later
        self.gaps=np.where(self.trace_diff>self.induction_gap)[0]
        print ("gaps in recordings", self.gaps)
        if len(self.gaps)==0:
            print ("whoops, no gap, possibly no induction")
        else:
            self.induction=self.gaps[-1]
            #if gap is position/array index n, then the gap is between n and n+1 in the tracetime array
            #column n is last trace of baseline, column n+1 contains first trace of follow-on
            #To analyze baseline, uses traces from n-(baseline_minutes*2) to n 
            #To analyze after induction, use traces from n+1 to the end
    
        if self.induction<self.pretraces:
            self.pretraces=self.induction #should this be -1 or line 131 +1 needed?
        self.starttrace=self.induction-self.pretraces+1 #first usable trace is pretraces prior to induction gap
        self.num_usable_traces=min(self.numtraces-self.starttrace,(two_hours+self.baseline_minutes)*self.tracesPerMinute)
        self.num_usable_minutes=int(np.ceil(float(self.num_usable_traces)/self.tracesPerMinute))
        self.goodtraces=range(self.starttrace,self.starttrace+self.num_usable_traces)
        print ("traces",self.numtraces,"usable", self.num_usable_traces,"start", self.starttrace,"two hours", (two_hours+self.baseline_minutes)*self.tracesPerMinute,"last usable", self.goodtraces[-1])
        #initialize arrays to hold measurements  
        self.base=np.zeros((self.num_usable_traces))
        self.peak=np.zeros((self.num_usable_traces))
        self.peaktime=np.zeros((self.num_usable_traces))
        self.pospeak=np.zeros((self.num_usable_traces))
        self.pospeaktime=np.zeros((self.num_usable_traces))
        self.FVsize=np.zeros((self.num_usable_traces))
        self.amp=np.zeros((self.num_usable_traces))
        self.baseline_start=np.zeros((self.num_usable_traces), dtype=int)
        self.popspikestart=np.zeros((self.num_usable_traces))
           
        return
        
    def findPopSpike(self):
        ######## Place for improvements.  
		## A. find artifact when it is small
		## B. find negative peak with less dependence on artifact decay time+FV width
		##     e.g., if artifact is positive, then could just find most negative point after artifact.
		##     will not work if artifact is negative
		##     alternative: low pass filter? would that eliminate artifact.
		## C. find fiber volley with less dependence on FV width
        for tracenum in self.goodtraces:
            index=tracenum-self.starttrace
            ######### Find artifact, earliest point where trace exceeds tempbaseline by some threshold
            artifactbegin=0
            tempbaseline=np.mean(self.Vm_traces[tracenum,0:self.baselinepoints])
            over_thresh=np.abs(self.Vm_traces[tracenum]-tempbaseline)>self.artifactthreshold
            if(np.all(~over_thresh)):
                artifactbegin = 0
            else:
                artifactbegin=np.min(np.where(over_thresh))
            self.baseline_start[index]=artifactbegin-int(sec_per_msec/self.dt) 
            artifact_end=int(artifactbegin+self.artifactdecay)
            ps_start=artifact_end+self.FVwidth_pts
            self.popspikestart[index]=ps_start*self.dt#artifact_end*self.dt
            #print(artifactbegin, self.latest_artifact, self.baseline_start[index], artifact_end, self.popspikestart[index])
            if artifactbegin==0 or artifactbegin>self.latest_artifact:
                print ("WARNING!!!!! NO ARTIFACT FOUND for trace",tracenum,'artifact', artifactbegin*self.dt, self.Vm_traces[tracenum,artifactbegin]-tempbaseline)
                self.base[index]=np.nan
                self.peaktime[index]=np.nan
                self.peak[index]=np.nan
                self.pospeaktime[index]=np.nan
                self.pospeak[index]=np.nan
                self.amp[index]=np.nan
                self.FVsize[index]=np.nan 
                #append data to problem traces dictionary for plotting later                   
                self.problemtraces.tracenum.append(tracenum)
                self.problemtraces.artifactbegin.append(artifactbegin)
                self.problemtraces.tempbaseline.append(tempbaseline)
            else:
                #Artifact found, find negative peak, baseline and the positive peak
                #Calculate base as mean of 1 ms prior to artifact
                self.base[index]=np.mean(self.Vm_traces[tracenum,self.baseline_start[index]:artifactbegin])
                #find negative peak and peak location, beginning at artifactdecay+FVwidth
                neg_peakpoint=self.Vm_traces[tracenum,ps_start:self.first_peak_end].argmin()+ps_start
                self.peaktime[index]=neg_peakpoint*self.dt
                #print('neg_peakpoint: ',neg_peakpoint,'tracenum: ',tracenum,ps_start,artifactbegin)
                self.peak[index]=self.Vm_traces[tracenum,neg_peakpoint]
                #find the peak which divides fiber volley and popspike
                #if the negative peak is found as part of artifact decay, plot as problem trace
                if (neg_peakpoint-ps_start)<2:
                    print ("WARNING!!!!! negative peak found at end of artifact decay + FVwidth for trace",tracenum,'artifact:', artifactbegin*self.dt, 'peak:', neg_peakpoint*self.dt)
                    self.pospeak[index]=np.nan
                    self.amp[index]=np.nan
                    self.FVsize[index]=np.nan
                    #append data to problem traces dictionary for plotting later
                    self.problemtraces.tracenum.append(tracenum)
                    self.problemtraces.artifactbegin.append(artifactbegin)
                    self.problemtraces.tempbaseline.append(tempbaseline)
                else:
                    pospeakpoint=self.Vm_traces[tracenum,artifact_end:neg_peakpoint].argmax()+artifact_end
                    self.pospeaktime[index]=pospeakpoint*self.dt
                    self.pospeak[index]=self.Vm_traces[tracenum,pospeakpoint]
                    #find fiber volley size - predefined width
                    if self.FVwidth_pts>0:
                        self.FVsize[index]=np.abs(self.Vm_traces[tracenum,artifact_end:ps_start].min()-self.base[index])
                    else:
                        self.FVsize[index]=0
                    #If pospeak is large noise artifact, this is problem trace
                    if (self.pospeak[index]-self.base[index])>self.noisethresh or (pospeakpoint-artifact_end)<2:
                        if (self.pospeak[index]-self.base[index])>self.noisethresh:
                            print ("WARNING!!!! positive peak is really big for trace",tracenum)
                        elif (pospeakpoint-artifact_end)<2:
                            print ("WARNING!!!!! positive peakpoint found at end of artifact decay for trace",tracenum,'artifact:', artifactbegin*self.dt, 'peak:', pospeakpoint*self.dt  , "baset", self.baseline_start[index])
                        self.pospeak[index]=np.nan
                        self.amp[index]=np.nan
                        self.FVsize[index]=np.nan
                          
                        self.problemtraces['tracenum'].append(tracenum)
                        self.problemtraces['artifactbegin'].append(artifactbegin)
                        self.problemtraces['tempbaseline'].append(tempbaseline)
                    else:
                        #calculate amplitude as difference between negative peak and either baseline or positive peak
                        if self.base[index] > self.pospeak[index]:
                            self.amp[index]=(self.base[index]-self.peak[index])
                        else: 
                            self.amp[index] = (self.pospeak[index]-self.peak[index])          
        return     
        
    def plotProblemTraces(self):
        spread = 0.1
        fig,axes=pyplot.subplots()
        fig.canvas.set_window_title('Problem Traces for '+params.exper)
        axes.set_ylabel('Vm (mV)')
        axes.set_xlabel('Time (sec), black | marks end of artifact decay, magenta diamond=')
        
        for tracenum,artifactbegin in zip(self.problemtraces["tracenum"],self.problemtraces["artifactbegin"]):
            index=tracenum-self.starttrace
            offset=0.001-self.baseline_start[index]*self.dt  #display 1msec prior to baseline_start
            artifact_end=int(artifactbegin+self.artifactdecay)
            axes.plot(self.time+offset,self.Vm_traces[tracenum,:]+index*0.1,label=tracenum)
            if not np.isnan(self.pospeaktime[index]):
                pospeakpoint = int(self.pospeaktime[index]/self.dt)
                axes.plot(self.pospeaktime[index]+offset,self.Vm_traces[tracenum,pospeakpoint]+index*spread,'mo')
            if not np.isnan(self.peaktime[index]):
                neg_peakpoint = int(self.peaktime[index]/self.dt)      
                axes.plot(self.peaktime[index]+offset,self.Vm_traces[tracenum,neg_peakpoint]+index*spread,'k*')
            axes.plot(self.popspikestart[index]-self.FVwidth+offset,self.Vm_traces[tracenum,artifact_end]+index*spread,'k|',ms=12)           
        axes.legend(fontsize=8, loc='best')           
        return 
    
    def plotGoodTraces(self):
        traces_per_panel = 40
        peak_decay = 0.016
        spread = 0.1
        panels=int(np.ceil(len(self.peaktime)/traces_per_panel))
        fig,axes=pyplot.subplots(1,panels)
        fig.canvas.set_window_title('Good traces for '+self.params.exper)
        for panel in range(0,panels):
            #axes[panel].clear()
            trace1=panel*traces_per_panel
            tracen=min(len(self.peaktime),(panel+1)*traces_per_panel) #what if the number of traces are under, fail,so u want to do this line
            for index in range(trace1,tracen):
                if ~np.isnan(self.base[index]): #some of these have nans in them, dont want to show those
                    st=self.baseline_start[index]
                    #st=0 ##### uncomment to display entire trace
                    end=int((self.peaktime[index]+peak_decay)/self.dt)
                    #end=len(Vm_traces[goodtraces[index]]) #####uncomment to display entire trace
                    offset=0.001-self.baseline_start[index]*self.dt  #display 1msec prior to baseline_start
                    axes[panel].plot(self.time[st:end]+offset,self.Vm_traces[self.goodtraces[index],st:end]+index*spread,label=self.goodtraces[index])
                    axes[panel].plot(self.peaktime[index]+offset,self.peak[index]+index*spread,'k*')
                    if np.isnan(self.pospeak[index]):
                         #print "nan detected", index, pospeaktime[index]+offset,int(pospeaktime[index]/dt)
                        axes[panel].plot(self.pospeaktime[index]+offset,self.Vm_traces[self.goodtraces[index],int(self.pospeaktime[index]/self.dt)]+index*spread,'mD')
                    elif self.base[index]>self.pospeak[index]:
                        axes[panel].plot(self.pospeaktime[index]+offset,self.pospeak[index]+index*spread,'r*')
                    else:
                        axes[panel].plot(self.pospeaktime[index]+offset,self.pospeak[index]+index*spread,'ro')
                    axes[panel].plot(self.dt*self.baseline_start[index]+offset,self.base[index]+index*spread,'bo')
                    axes[panel].plot(self.popspikestart[index]+offset,self.Vm_traces[self.goodtraces[index],int(self.popspikestart[index]/self.dt)]+index*spread,'b|')
                    axes[panel].plot(self.popspikestart[index]-self.FVwidth+offset,self.Vm_traces[self.goodtraces[index],int(self.popspikestart[index]/self.dt)]+index*spread,'k|')
            axes[panel].set_xlabel('Time (sec)')
            axes[panel].legend(fontsize=7, loc='right')
        fig.suptitle('black * =popspike , red o,* =pospeak, blue o =baseline, | =ps start/FV start')
        axes[0].set_ylabel('Vm (mV)')
        fig.canvas.draw()        
        return    
    
    def normalizeAmp(self):
        pretrace_amp=self.amp[0:self.pretraces]
        self.meanpre=pretrace_amp[~np.isnan(pretrace_amp)].mean()
        #meanpre=self.amp[~np.isnan(self.amp[0:pretraces])].mean()
        self.pretraces=self.baseline_minutes*self.tracesPerMinute
        #plot traces with extracted peaks to verify
        #Next, normalize amplitude by mean popspike from preinduction traces, and
        # average two samples per minute to produce a single value per minute, both for popspike and fibre volley
        self.popspikenorm=np.zeros((self.num_usable_minutes))
        self.FVnorm=np.zeros((self.num_usable_minutes))
        self.popspikeminutes=np.zeros((self.num_usable_minutes))
        self.meanFVpre=self.FVsize[0:self.pretraces].mean()
        for sample in range(self.num_usable_minutes):
            self.popspikenorm[sample]=np.nanmean(self.amp[sample*self.tracesPerMinute:(sample+1)*self.tracesPerMinute])/self.meanpre	
            self.popspikeminutes[sample]=(self.tracetime[self.starttrace+sample*self.tracesPerMinute]-self.tracetime[self.induction])/60	#60 sec per minute
            if self.FVwidth > 0:
                self.FVnorm[sample]=self.FVsize[sample*self.tracesPerMinute:(sample+1)*self.tracesPerMinute].mean()/self.meanFVpre
            else:
                self.FVnorm[sample]=0       
        return
    
    def summaryMeasure(self):
        #Extract a few measures for statistical analysis, e.g. extract the 30, 60, 90, 120 min PopSpike means and FiberVolley
        #Times are specified in sample_times array 
        self.popspike_timesamples=np.zeros((len(self.sample_times)+1))
        self.FV_timesamples=np.zeros((len(self.sample_times)+1))
        self.popspike_timesamples[0]=self.meanpre
        self.FV_timesamples[0]=self.meanFVpre
        deletetime=[]
        for index in range(0,len(self.sample_times)):
            norm_index=self.sample_times[index]+self.baseline_minutes
            if len(self.popspikenorm)<norm_index:
                print ("calculating time samples; recording of ", self.popspikeminutes[-1], "min shorter than ", self.sample_times[index])
                self.popspike_timesamples[index+1]=np.nan
                self.FV_timesamples[index+1]=np.nan
                deletetime.append(index)
            else:
                self.popspike_timesamples[index+1]=np.nanmean(self.popspikenorm[norm_index-self.sample_window:norm_index])
                self.FV_timesamples[index+1]=self.FVnorm[norm_index-self.sample_window:norm_index].mean()
        #May need to delete this line if column_stack in GrpAvgPopSpike complains about different length arrays
        self.sample_times=np.int_(np.delete(self.sample_times,deletetime)+(self.baseline_minutes-self.sample_window/2))        
        return
    
    def showSummaryPlot(self):
        print ("ready for summary plot, and fit line to baseline")
        #Plot popspikenorm vs minutes, fit line to baseline, print slope
        validbasepts=~np.isnan(self.popspikenorm) #[0:baseline_minutes])
        validbasepts[self.baseline_minutes:]=False           
        popt,pcov=optimize.curve_fit(psu.line,self.popspikeminutes[validbasepts],self.popspikenorm[validbasepts])
        Aopt,self.Bopt=popt
        Astd,self.Bstd=np.sqrt(np.diag(pcov))
        if  np.abs(self.Bopt)-self.slope_std_factor*self.Bstd>0:
            print( "********************** WARNING ************** \n Baseline Slope of ", self.Bopt,"+/- 2*",self.Bstd, "sig diff than zero")
        else:
            print( "**********************Baseline Seems to Be Fine********************: slope=", self.Bopt)
        fig,axes=pyplot.subplots()
        fig.canvas.set_window_title('Summary ')
        axes.plot(self.popspikeminutes,self.popspikenorm,'b.')
        axes.plot(self.popspikeminutes[0:self.baseline_minutes],psu.line(self.popspikeminutes[0:self.baseline_minutes],Aopt,self.Bopt),'r')
        start=np.insert(self.sample_times,0,0)
        time_summary=self.popspikeminutes[start]
        print ("summary at", len(start), "time points \n",np.column_stack((time_summary,self.popspike_timesamples[0:len(start)])))
        axes.plot(time_summary,self.popspike_timesamples[0:len(start)],'ko')
        axes.set_xlabel('Time (minutes)')
        axes.set_ylabel ('Normalized Pop Spike ')
        return        

    def saveData(self): 
        #tracedict includes all traces that you potentially want to analyze by groups in GroupAvg
        #modify argsdict to include all the other bits of information you want, and change dictionary keys to be meaningful for you.
        tracedict={'amp':self.amp,'peaktime':self.peaktime,'FVsize':self.FVsize,'popspikenorm':self.popspikenorm,'popspikeminutes':self.popspikeminutes, 
                   'FVnorm': self.FVnorm,'PS_mean':self.popspike_timesamples,'FV_means':self.FV_timesamples,'slope':self.Bopt,'slope_std':self.Bstd}
        #This writes the file for GrpAvg.py
        #it includes all the experiment parameters
        datadict = dict(trace=tracedict, #datadict name wont be accessible outside this program but keys are
                parameters=params, anal_params=self.anal_params) #first is key second is value
        outpickle = self.outputdir+self.params.exper + '.pickle'
        with open(outpickle, 'wb') as f:
            pickle.dump(datadict, f) #saves the contents of datadict in file "f"...       
        return
            
######## Main ########## ARGS here are for testing the code            
ARGS="20160822orangedm105 M 30 none 10.5 DM -decay 1.6 -FVwidth 1.1"            
try:
    commandline = ARGS.split() # ARGS="2015_04_07_02 M 30 none 10.5 DM "
    do_exit = False
except NameError: 
    commandline = sys.argv[1:]
    do_exit = True    
    
params=psu.parse_args(commandline,do_exit,0)

pop_spike = PopSpike(params) 
pop_spike.read_datafile()
pop_spike.find_induction_gap()
pop_spike.findPopSpike()
pop_spike.plotProblemTraces()
pop_spike.plotGoodTraces()
pop_spike.normalizeAmp()
pop_spike.summaryMeasure()
pop_spike.showSummaryPlot()
pop_spike.saveData()

#####