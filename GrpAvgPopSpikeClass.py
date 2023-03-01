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
        self.baseline_min=0.4
        self.baseline_max=2
        self.plot_individuals=int(params.plot_ctrl[1])  #to plot PopSpike vs time for each experiment in a group
        self.plot_corr=int(params.plot_ctrl[2]) #to plot  correlation between LTP (at summarytime) and age or baseline epsp
        self.print_info=1
        self.colors=pyplot.get_cmap('viridis')
        self.color_incr=32
        self.BAD = []
       
    def read_data(self): 
        DATAS = []
        PARAMS = pd.DataFrame()
        ANAL = pd.DataFrame()
        for outfname in self.outfnames:
            with open(outfname, 'rb') as f:
                datadict = pickle.load(f)
            if self.print_info and 'slope' in datadict['trace'].keys():
                print ("file read:", datadict['parameters'],", baseline slope=", round(datadict['trace']['slope'],6))
            exper_param = datadict['parameters']
            numnan=sum(np.isnan(datadict['trace']['popspikenorm']))
            ############ Select subset of files based on user specified criteria ##########
            drug_parts=exper_param.drug.split('_')
            drug_name=drug_parts[0]
            if len(drug_parts)>1:
                drug_conc=drug_parts[1]
            ignore1 = ((self.params.sex and self.params.sex != exper_param.sex) or 
                          (self.params.age is not None and self.params.age >= exper_param.age) or
                          (self.params.maxage is not None and self.params.maxage <= exper_param.age) or
                          (self.params.drug and self.params.drug != drug_name) or
                          (self.params.region and self.params.region != exper_param.region) or
                          (self.params.theta and self.params.theta != exper_param.theta))
            ########## identify experiments that do not meet includsion criteria ##########
            exclude = ((len(datadict['trace']['popspikenorm'])<self.minimum_sweeps) or
                          (np.abs(datadict['trace']['slope'])>self.slope_threshold) or
                          (np.abs(datadict['trace']['slope'])-self.slope_std_factor*datadict['trace']['slope_std']>0) or
                          numnan>self.nan_threshold or #do not use trace if too many missing popspikes
						  datadict['trace']['PS_mean'][0]<self.baseline_min)
            if exclude and not ignore1:
                if len(datadict['trace']['popspikenorm'])<self.minimum_sweeps:
                    print ("!!!NOT ENOUGH goodtraces", datadict['parameters'].exper, len(datadict['trace']['popspikenorm']))
                elif np.abs(datadict['trace']['slope'])>self.slope_threshold or (np.abs(datadict['trace']['slope'])-self.slope_std_factor*datadict['trace']['slope_std']>0):
                    print ("!!!BAD baseline", datadict['parameters'].exper,  "slope u,s", round(datadict['trace']['slope'],6),round(datadict['trace']['slope_std'],6))
                elif datadict['trace']['PS_mean'][0]<self.baseline_min:
                    print ("!!!PSP amp too low", datadict['trace']['PS_mean'][0])
                else:
                    print('!!! Too many Nans',numnan)
                self.BAD.append(datadict) #Don't change until you have the ploting function
                #text=raw_input('continue? (y/n)')
                #Ignore 2 is T if pop spike 2 doesn't have enough sweep. Ignore 2 means do not use it if its bad 
                if self.print_info:
                    print ("***** excluding:", datadict['parameters'])
                next 
            elif ignore1: #regardless of whether data meets exclusion criteria
                if self.print_info:
                    print ("***** ignoring:", datadict['parameters'])
                next 
            else: #data meets all criteria
                if np.isnan(datadict['trace']['popspikenorm']).any(): #inform about nans, even if using data
                    print ("########## np.nan detected", numnan,'times in',exper_param.exper, datadict['anal_params'])
                if datadict['trace']['popspikenorm'][-1]==0.0: #last popspike is 0, why????
                    print ("@@@@@@@@@ check popSpikeAnal for", exper_param.exper, datadict['anal_params'])
                DATAS.append(datadict['trace'])
                params=vars(datadict['parameters'])
                params['slope_std_factor']=params['slope_std']
                del params['slope_std']
                params['date']=params['exper'][0:8]
                PARAMS = PARAMS.append(params, ignore_index = True)
                ANAL=ANAL.append(datadict['anal_params'],ignore_index=True)
        
        df2 = pd.DataFrame(DATAS, dtype=float)
        self.whole_df=pd.concat([PARAMS,ANAL, df2], axis = 1)
        #### combine experiments if same conditions collected from same mouse
        print('####################################### ') 
        avg_col=['popspikeminutes','popspikenorm','PS_mean','FV_means','FVnorm','FVsize','peaktime','slope_std','slope']
        grp_params=['date','region','theta','sex','drug','age']
        params=	grp_params+['exper','criticaltimes']#vars(datadict['parameters']).keys()
        grp_data=self.whole_df.groupby(grp_params)
        all_new_data=[];drop_indices=[]
        for grp in grp_data.groups.keys():
             if grp_data.get_group(grp).age.count()>1:
                  new_data={}
                  drop_indices.append(grp_data.get_group(grp).index)
                  for ky in avg_col:
                     sizes=[self.whole_df[ky].iloc[x].size for x in drop_indices[-1]]
                     size_std=np.std(sizes)
                     if  size_std==0:
                          new_data[ky]=grp_data.get_group(grp)[ky].mean() #replace with nanmean somehow?
                     else:
                          new_data[ky],_,_=grp_utl.exp_avg(grp_data.get_group(grp)[ky]) #replace with nanmean somehow?
                  for ky in params:
                      new_data[ky]=grp_data.get_group(grp)[ky].values[0]
                  all_new_data.append(new_data) 
        #print('before dropping,', drop_indices, 'whole_df length=',len(self.whole_df),"\n",self.whole_df.exper)	  
        for ind in [dp for dpset in drop_indices for dp in dpset]:
            self.whole_df=self.whole_df.drop(ind) 
            print('drop_index',ind,'new df length',len(self.whole_df))
        for new_data in all_new_data:
            self.whole_df = self.whole_df.append(new_data, ignore_index = True)
        print (len(self.whole_df),'experiments met criteria.')
    
    def group_data(self, sepvarlist, newcolumnname=None): 
        if newcolumnname is not None:
            self.whole_df[newcolumnname]=self.whole_df.apply(lambda row: row.drug.split('_')[0],axis=1)
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
            self.filenm[grp]=self.construct_filename(sepvarlist,grp)
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
        axes[0].set_ylabel('popspike-bad slope or psp amp')
        axes[1].set_ylabel('popspike-nan')
        axes[1].set_xlabel('Time (min)')
        for bad_data in self.BAD: 
            if np.isnan(bad_data['trace']['popspikenorm']).any():
             	  ### instead of plotting here, should put in separate container to plot later
                  print ("########## np.nan detected", bad_data['parameters'].exper, bad_data['anal_params'])
                  p=axes[1].plot(bad_data['trace']['popspikeminutes'],bad_data['trace']['popspikenorm'],'+', label=bad_data['parameters'].exper)
                  color = p[-1].get_color()
                  nan_index=np.argwhere(np.isnan(bad_data['trace']['popspikenorm']))
                  axes[1].plot(bad_data['trace']['popspikeminutes'][nan_index],np.ones((len(nan_index))),'o', color=color)
            else:
                  axes[0].plot(bad_data['trace']['popspikeminutes'],bad_data['trace']['popspikenorm'],'.', label=bad_data['parameters'].exper)
            if bad_data['trace']['popspikenorm'][-1]==0.0:
                    print ("@@@@@@@@@ check popSpikeAnal for", bad_data['parameters'].exper, bad_data['anal_params'])
                #print 'OK: {}'.format(exper_param)
            if bad_data['trace']['popspikeminutes'][1]-bad_data['trace']['popspikeminutes'][0]<0.8:
                    print(bad_data['parameters'].exper,bad_data['anal_params']['artdecay'], bad_data['anal_params']['FVwidth'],bad_data['trace']['popspikeminutes'])       
        axes[0].legend(fontsize=8, loc='best')
        axes[1].legend(fontsize=8, loc='best')
        fig.canvas.draw()

    def write_stat_data(self, newcolumnname):
        params = ['exper', 'sex', 'age', 'region', 'theta', 'drug', newcolumnname]
        SASoutput = self.whole_df[params] #in to a 2d array you write it into SASoutput. 
        SASheader= ' '.join(params)
        sample_times=[30,60] ####### FIXME: hard coded sample times
        for row in range(len(self.whole_df)):
            for index in range(0,len(sample_times)):
                norm_index=np.min(np.where(self.whole_df.iloc[row]['popspikeminutes']>sample_times[index]))
                self.whole_df.iloc[row].PS_mean[index+1]=np.nanmean(self.whole_df.iloc[row]['popspikenorm'][norm_index-2:norm_index+1])
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
    def bar_graph_data(self,exclude_name):
        lines=[]
        for grp in Grp_PS.grp_data.groups.keys():
           #[0:3] gives you baseline amp, 30min and 60 min
           PS_mean=Grp_PS.grp_data.get_group(grp).PS_mean.values
           line=[grp]+[np.mean(PS_mean,axis=0)[1]*100]+[(np.std(PS_mean,axis=0)[1]/np.sqrt(len(PS_mean)))*100]
           lines.append(line)
           #print(line,np.std(PS_mean,axis=0)[1])
        prefix=self.common_filename()
        f=open(prefix+"_BarGraphmeans.txt", 'w')
        f.write('  '.join(Grp_PS.sepvarlist)+'mean    sterr \n')
        np.savetxt(f, lines, fmt='%s', delimiter='   ')
        f.close()
        for grp in Grp_PS.grp_data.groups.keys():
            PS_mean=Grp_PS.grp_data.get_group(grp).PS_mean.values
            data_column=[ps[1]*100 for ps in PS_mean] #convert the 30 min time point to percent
            f=open(self.filenm[grp]+"_points.txt",'w')
            if isinstance(grp,tuple):
                columnheading=[prefix+str(p) for i,p in enumerate(grp) if self.sepvarlist[i] not in exclude_name]     #assumes you don't want theta as part of column name ; assumes drug comes before sex in ARGS           
            else:
                columnheading=[prefix+grp]
            f.write('_'.join(columnheading)+"_data\n")
            np.savetxt(f,data_column, fmt='%7.5f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
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
    def common_filename(self):			
        filnm=''
        for p in [self.params.drug,self.params.sex,self.params.region,self.params.theta]:
            if p is not None:
                filnm=filnm+str(p)+'_'
        return filnm
    def construct_filename(self,sepvarlist,paramgrp): 
        filnm=self.common_filename()
        if len(sepvarlist)==1:
           paramgrp=[paramgrp]
        filnm=filnm+'_'.join([sv+str(v) for sv,v in zip(sepvarlist,paramgrp)])
        return filnm
        
#ARGS = "-outputdir ../pickle/ -sepvarlist sex region theta -plot_ctrl 000"      

if __name__ =='__main__':        
	exclude_name=['theta'] #use to exclude variable(s) from column name in _points files
	        
	try:
	 	commandline = ARGS.split() #in python: define space-separated ARGS string
	 	do_exit = False
	except NameError: #if you are not in python, read in filename and other parameters
	 	commandline = sys.argv[1:]
	 	do_exit = True
	
	params=psu.parse_args(commandline,do_exit,1)
	#print('params=',params,'commandline =',commandline)            
	Grp_PS=GrpPopSpike(params)
	Grp_PS.pattern = Grp_PS.subdir+'*.pickle'
	Grp_PS.outfnames = sorted(glob.glob(Grp_PS.pattern))
	Grp_PS.outfnames = [f for f in Grp_PS.outfnames if not f.endswith('IO.pickle')]
	if not len(Grp_PS.outfnames):
	    sys.exit('************ No files found *********')
	if len(Grp_PS.outfnames):
	    Grp_PS.read_data()
	    Grp_PS.plot_bad()
	    newcolumnname='drug_combined'
	    Grp_PS.group_data(Grp_PS.sepvarlist,newcolumnname) 
	    if int(params.plot_ctrl[0])>0:
	        plot_cols=int(params.plot_ctrl[0])
	    else:
	        plot_cols=None
	    grp_utl.plot_groups(Grp_PS.avg_PS,Grp_PS.stderr_PS,Grp_PS.minutes,Grp_PS.samples,Grp_PS.filenm,Grp_PS.sepvarlist,Grp_PS.sepvardict,plot_cols)
	    Grp_PS.write_traces() 
	    Grp_PS.write_stat_data(newcolumnname) 
	    Grp_PS.bar_graph_data(exclude_name)
	##### After working, deal with continuous valued groups, e.g. age and light level
	### update GrpAvgPSP
	

              
                
                
                
                
                
                
                
                