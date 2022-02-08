# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:34:30 2021

@author: nehas
"""
import numpy as np
import sys  
import pandas as pd 
import glob
import pickle
import os
print (os.getcwd())
from matplotlib import pyplot
import pop_spike_utilities as psu
import GrpPlotUtil2 as grp_utl
from scipy import optimize

class GrpPopSpike:
    def __init__(self, params): 
        self.params = params
        self.subdir = params.outputdir
        self.summarytime = params.samp_time
        self.slope_std_factor=params.slope_std 
        self.sepvarlist= params.sepvarlist
        #additional parameters.  consider making commandline parameters
        self.minimum_sweeps=45     
        self.slope_threshold=0.01  
        self.nan_threshold=20 
		self.baseline_amp=0.5
        self.plot_individuals=int(params.plot_ctrl[1])  #to plot PopSpike vs time for each experiment in a group
        self.plot_corr=int(params.plot_ctrl[2]) #to plot  correlation between LTP (at summarytime) and age or baseline epsp
        self.print_info=1
        self.colors=pyplot.get_cmap('viridis')
        self.color_incr=32
        self.pattern = self.subdir+'20*.pickle'
        self.outfnames = sorted(glob.glob(self.pattern))
        if not len(self.outfnames):
            print('************ No files found *********')
        self.BAD = []
       
    def read_data(self): 
        DATAS = []
        PARAMS = pd.DataFrame()
        ANAL = pd.DataFrame()
        for outfname in self.outfnames:
            with open(outfname, 'rb') as f:
                datadict = pickle.load(f)
            if self.print_info:
                print ("file read:", datadict['parameters'],", baseline slope=", round(datadict['trace']['slope'],6))
            exper_param = datadict['parameters']
            numnan=sum(np.isnan(datadict['trace']['popspikenorm']))
            ############ Select subset of files based on user specified criteria ##########
            ignore1 = ((self.params.sex and self.params.sex != exper_param.sex) or 
                          (self.params.age is not None and self.params.age >= exper_param.age) or
                          (self.params.maxage is not None and self.params.maxage <= exper_param.age) or
                          (self.params.drug and self.params.drug != exper_param.drug) or
                          (self.params.region and self.params.region != exper_param.region) or
                          (self.params.theta and self.params.theta != exper_param.theta))
            ########## identify experiments that do not meet includsion criteria ##########
            ignore2 = ((len(datadict['trace']['popspikenorm'])<self.minimum_sweeps) or
                          (np.abs(datadict['trace']['slope'])>self.slope_threshold) or
                          #-self.slope_std_factor*datadict['trace']['slope_std']
                          numnan>self.nan_threshold) or
						  datadict['trace']['popspike_timesamples'][0]<self.baseline_amp) #do not use trace if too many missing popspikes
            if ignore2 and not ignore1:
                if len(datadict['trace']['popspikenorm'])<self.minimum_sweeps:
                    print ("!!!NOT ENOUGH goodtraces", datadict['parameters'].exper, len(datadict['trace']['popspikenorm']))
                elif np.abs(datadict['trace']['slope'])>self.slope_threshold:
                    print ("!!!BAD baseline", datadict['parameters'].exper,  "slope u,s", round(datadict['trace']['slope'],6),round(datadict['trace']['slope_std'],6))
                else:
                    print('!!! Too many Nans',numnan)
                self.BAD.append(datadict) #Don't change until you have the ploting function
                #text=raw_input('continue? (y/n)')
                #Ignore 2 is T if pop spike 2 doesn't have enough sweep. Ignore 2 means do not use it if its bad 
            if ignore1 or ignore2:
                if self.print_info:
                    print ("***** ignoring:", ignore1,ignore2,datadict['parameters'])
                next 
            else:
                if np.isnan(datadict['trace']['popspikenorm']).any():
                    print ("########## np.nan detected", numnan,'times in',exper_param.exper, datadict['anal_params'])
                if datadict['trace']['popspikenorm'][-1]==0.0:
                    print ("@@@@@@@@@ check popSpikeAnal for", exper_param.exper, datadict['anal_params'])
                    
                    ###### change this to put into PANDAS dataframe
            # Wait until you can accumilate the data. 
            #Take out the selfs in datadict b/c it changes.
                DATAS.append(datadict['trace'])
                PARAMS = PARAMS.append(vars(datadict['parameters']), ignore_index = True)
                ANAL=ANAL.append(datadict['anal_params'],ignore_index=True)
        
        df2 = pd.DataFrame(DATAS, dtype=float)
        self.whole_df=pd.concat([PARAMS,ANAL, df2], axis = 1)
        print (len(self.whole_df),'experiments met criteria.')
    
    def group_data(self, sepvarlist): 
        self.grp_data = self.whole_df.groupby(sepvarlist)
        self.avg_PS={}
        self.stderr_PS={}
        self.samples={}
        self.minutes={}
        self.filenm={}
        for grp in self.grp_data.groups.keys(): #consider changing tuple to string using join
            avg,std,count=grp_utl.exp_avg(self.grp_data.get_group(grp).popspikenorm)
            self.avg_PS[grp]=avg
            self.stderr_PS[grp]=std/np.sqrt(len(self.grp_data.get_group(grp)))
            self.samples[grp]=count
            self.minutes[grp]=grp_utl.exp_avg(self.grp_data.get_group(grp).popspikeminutes)[0]
            #these don't work because popspikenorm not always the same length
            #self.avgpopspikenorm[grp]=np.nanmean(self.grp_data.get_group(grp).popspikenorm) 
            #self.minutes[grp]= np.mean(self.grp_data.get_group(grp).popspikeminutes) 
            #self.stderrpopspikenorm[grp]= np.std(self.grp_data.get_group(grp).popspikenorm)/np.sqrt(self.samples(grp)) 
            self.filenm[grp]=grp_utl.construct_filename(sepvarlist,grp)
        #this dictionary is used for plotting in GrpPlotUtil
        self.sepvardict={}
        for sv in self.sepvarlist:
            self.sepvardict[sv]=np.unique(self.whole_df[sv])
  
    def plot_bad(self):
        if len(self.BAD):
            for bad_data in self.BAD:
                if bad_data['trace']['slope']>self.slope_threshold:
                    slopes=[];slopestd=[]
                    for stpt in [0,5]:
                        baseline_end=np.min(np.where(bad_data['trace']['popspikeminutes']>0))
                        validbasepts=np.where(~np.isnan(bad_data['trace']['popspikenorm'][stpt:baseline_end]))[0]+stpt
                        popt,pcov=optimize.curve_fit(psu.line,bad_data['trace']['popspikeminutes'][validbasepts],bad_data['trace']['popspikenorm'][validbasepts])
                        slopes.append(popt[1])
                        slopestd.append(np.sqrt(np.diag(pcov))[1])
                    print ("bad baseline", bad_data['parameters'].exper, "slope",round(slopes[0],5), "+/-", round(slopestd[0],5), "10min slope", round(slopes[1],5), round(slopestd[1],5))

        fig,axes=pyplot.subplots(2,1)
        fig.canvas.set_window_title('Problems')
        axes[0].set_ylabel('popspike-bad slope')
        axes[1].set_ylabel('popspike-nan')
        axes[1].set_xlabel('Time (min)')
        for bad_data in self.BAD: 
            axes[0].plot(bad_data['trace']['popspikeminutes'],bad_data['trace']['popspikenorm'],'.', label=bad_data['parameters'].exper)
            if np.isnan(bad_data['trace']['popspikenorm']).any():
             	  ### instead of plotting here, should put in separate container to plot later
                  print ("########## np.nan detected", bad_data['parameters'].exper, bad_data['anal_params'])
                  axes[1].plot(bad_data['trace']['popspikeminutes'],bad_data['trace']['popspikenorm'],'+', label=bad_data['parameters'].exper)
                  nan_index=np.argwhere(np.isnan(bad_data['trace']['popspikenorm']))
                  axes[1].plot(bad_data['trace']['popspikeminutes'][nan_index],np.ones((len(nan_index))),'o', label= bad_data['parameters'].exper)
            if bad_data['trace']['popspikenorm'][-1]==0.0:
                    print ("@@@@@@@@@ check popSpikeAnal for", bad_data['parameters'].exper, bad_data['anal_params'])
                #print 'OK: {}'.format(exper_param)
            if bad_data['trace']['popspikeminutes'][1]-bad_data['trace']['popspikeminutes'][0]<0.8:
                    print(bad_data['parameters'].exper,bad_data['anal_params']['artdecay'], bad_data['anal_params']['FVwidth'],bad_data['trace']['popspikeminutes'])       
        axes[0].legend(fontsize=8, loc='best')
        axes[1].legend(fontsize=8, loc='best')
        fig.canvas.draw()

    def write_stat_data(self):
        params = ['exper', 'sex', 'age', 'region', 'estradiol', 'theta', 'drug']
        SASoutput = self.whole_df[params] #in to a 2d array you write it into SASoutput. 
        SASheader= ' '.join(params)
        for col in range(len(self.whole_df.PS_mean[0])):
            psmean=[row[col] for row in self.whole_df.PS_mean]
            if not np.all([np.isnan(k) for k in psmean]):
                SASoutput=np.column_stack((SASoutput,psmean))
                if col==0:
                    SASheader += ' baseline'
                else:
                    SASheader += ' normpopspike'+str(col*30)
        f=open("PARAMSforSAS.txt", 'w')
        f.write(SASheader +"\n")
        np.savetxt(f, SASoutput, fmt='%s', delimiter='   ')
        f.close()

    def write_traces(self):
        frac_2_percent = 100
        for gnum in self.grp_data.groups.keys():
            f=open(self.filenm[gnum]+".txt",'w')                      
            outputdata=np.column_stack((self.minutes[gnum], self.samples[gnum],
                                frac_2_percent*self.avg_PS[gnum], frac_2_percent*self.stderr_PS[gnum]))
            header="time "+self.filenm[gnum]+"count "+self.filenm[gnum]+"normpsAVG "+self.filenm[gnum]+"normpsSE "
            f.write(header+"\n")
            np.savetxt(f,outputdata, fmt='%7.5f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
            f.close()
        
ARGS = "-outputdir ../pickle/ -sepvarlist sex region theta -plot_ctrl 000"      
        
try:
 	commandline = ARGS.split() #in python: define space-separated ARGS string
 	do_exit = False
except NameError: #if you are not in python, read in filename and other parameters
 	commandline = sys.argv[1:]
 	do_exit = True

params=psu.parse_args(commandline,do_exit,1)
print(params, commandline)            
Grp_PS=GrpPopSpike(params)
if len(Grp_PS.outfnames):
    Grp_PS.read_data()
    Grp_PS.plot_bad()
    Grp_PS.group_data(Grp_PS.sepvarlist) 
    if int(params.plot_ctrl[0])>0:
        plot_cols=int(params.plot_ctrl[0])
    else:
        plot_cols=None
    grp_utl.plot_groups(Grp_PS.avg_PS,Grp_PS.stderr_PS,Grp_PS.minutes,Grp_PS.samples,Grp_PS.filenm,Grp_PS.sepvarlist,Grp_PS.sepvardict,plot_cols)
    Grp_PS.write_traces() 
    Grp_PS.write_stat_data() 

##### After working, deal with continuous valued groups, e.g. age and light level
### update GrpAvgPSP
  

              
                
                
                
                
                
                
                
                