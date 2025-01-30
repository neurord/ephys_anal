# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:34:30 2021

@author: nehas
"""
import numpy as np
import sys  
import pandas as pd 
import glob
import os
print (os.getcwd())
from matplotlib import pyplot
import ArgParser as argp
import PatchAnalPlots as pu
import GrpPlotUtil2 as grp_utl

SEC_PER_MIN=60 #global

class GrpPatch:
    def __init__(self, params): 
        self.params = params
        self.subdir = params.outputdir
        self.slope_std_factor=params.slope_std 
        self.sepvarlist= params.sepvarlist
        self.sample_times=params.samp_time
        self.window=params.window
        self.plot_individuals=int(params.plot_ctrl[1])  #to plot PSP vs time for each experiment in a group
        self.plot_corr=int(params.plot_ctrl[2]) #to plot correlation between LTP and age or baseline epsp
        #additional parameters.  FIXME: add to arg parser
        self.minimum_sweeps=10#20 #5 min pre and 15 min follow-up
        self.slope_threshold=0.01
        self.nan_threshold=10 
        self.baseline_min=0#0.4
        self.baseline_max=2
        self.print_info=1 
        self.single_params=['region','genotype','age','drug','ID','time_to_induct','pre_num']
        self.fname_vars=['region','genotype','sex','drug']
        self.IVIF_variables=['Im','Vm','latency','num_spikes'] #also, APheight, APwidth, AHP_amp, AHP_time

    def read_data(self): 
        self.BAD = {}
        DATAS = []
        PARAMS = []
        ANAL_PARAMS=[]
        IVIFset=[]
        for outfname in self.outfnames:
            with open(outfname, 'rb') as f:
                datadict = np.load(f,allow_pickle=True)
                data=datadict['data'].item()
                exper_param=datadict['params'].item()
                traces=datadict['trace'].item()
                IV_IF=datadict['IV_IF'].item()
                anal_params=datadict['anal_params'].item() #use if need to re-analyze experiment
                celltype=exper_param['celltype']
                print_params={k:v for k,v in exper_param.items() if not isinstance(v,dict) or k=='celltype'}
            if self.print_info and 'slope' in data.keys():
                print ("file read:", print_params,", baseline slope=", [round(sl,6) for sl in data['slope'].values()])
            ############ Select subset of files based on user specified criteria ##########
            ignore = ((self.params.sex and self.params.sex != exper_param['sex']) or 
                          (self.params.age is not None and self.params.age >= exper_param['age']) or
                          #(self.params.maxage is not None and self.params.maxage <= exper_param['age']) or
                          (self.params.drug and self.params.drug != exper_param['drug']) or
                          (self.params.region and self.params.region != exper_param['region']))# or
                          #(self.params.induction and self.params.induction !=exper_param['induction']['induct_freq']))
            ########## identify experiments that do not meet includsion criteria ##########
            for hs in celltype.keys():
                numnan=sum(np.isnan(traces['normPSP'][hs]))
                exclude = ((len(traces['normPSP'][hs])<self.minimum_sweeps) or
                                (np.abs(data['slope'][hs])>self.slope_threshold) or
                                (np.abs(data['slope'][hs])-self.slope_std_factor*data['slope_std'][hs]>0) or
                                numnan>self.nan_threshold or #do not use trace if too many missing popspikes
                                np.mean(data['meanpre'][hs])<self.baseline_min) #do not use if baseline psp is too small
                if exclude and not ignore:
                    if len(traces['normPSP'][hs])<self.minimum_sweeps:
                        print ("!!!NOT ENOUGH goodtraces", exper_param['exper'], hs, len(traces['normPSP'][hs]))
                    elif np.abs(data['slope'][hs])>self.slope_threshold or (np.abs(data['slope'][hs])-self.slope_std_factor*data['slope_std'][hs]>0):
                        print ("!!!BAD baseline", exper_param['exper'],  hs, "slope u,s", round(data['slope'][hs],6),round(data['slope_std'][hs],6))
                    elif data['meanpre'][hs]<self.baseline_min:
                        print ("!!!PSP amp too low", data['meanpre'][hs])
                    else:
                        print('!!! Too many Nans',exper_param['exper'], hs, numnan)
                    self.BAD[exper_param['exper']+'_'+hs]={'normPSP':traces['normPSP'][hs],'psptime':traces['psptime'][hs],'slope':data['slope'][hs],'slope_std':data['slope_std'][hs],
                                                        'params':{'celltype':celltype[hs],'age':exper_param['age'],'drug':exper_param['drug'],'region':exper_param['region']}}
                    #Exclude is True if not enough sweeps or otherwise bad data.  Ignore is true if not the subset of data desired
                    if self.print_info:
                        print ("***** excluding:", exper_param['exper']+'_'+hs)
                    next 
                elif ignore: #regardless of whether data meets exclusion criteria
                    if self.print_info:
                        print ("***** ignoring:", exper_param['exper']+'_'+hs)
                    next 
                else: #data meets all criteria
                    if np.isnan(traces['normPSP'][hs]).any(): #inform about nans, even if using data
                        print ("########## np.nan detected", numnan,'times in',exper_param['exper']+'_'+hs)
                    if traces['normPSP'][hs][-1]==0.0: #last popspike is 0, why????
                        print ("@@@@@@@@@ check PatchAnal for", exper_param['exper']+'_'+hs)
                    cell_param=self.extract_params(exper_param,hs,data)
                    PARAMS.append(cell_param)
                    celltrace={k:val[hs] for k,val in traces.items()}
                    celltrace['PSPsamples']=self.summary_measure(celltrace,cell_param)
                    DATAS.append(celltrace)                    
                    IVIFset.append(self.IVIFdata(IV_IF,hs,cell_param))
        
        dfdata = pd.DataFrame(DATAS)
        dfparams=pd.DataFrame(PARAMS)
        df_ivif=pd.DataFrame(IVIFset)
        self.whole_df=pd.concat([dfparams,dfdata,df_ivif], axis = 1)
        self.single_params.append('celltype')
        self.single_params.remove('pre_num')

    def IVIFdata(self,IVIF,hs,params):
        IVIF_data={}
        IFneeded=True
        for key,vals in IVIF['IV'][hs].items():
            #identify IV and IF from prior to baseline 
            routine=key.split('_')[0]
            if int(routine[1:])<params['pre_num']: #FIXME: can also get data after induction - different variable name
                if 'IV' in key:
                    for arr in self.IVIF_variables:
                        IVIF_data[arr]=vals[arr] #this assumes that IV is first
                if 'IF' in key and IFneeded:
                    IFneeded=False #prevent multiple IF curves from pre-induction
                    for arr in self.IVIF_variables:
                        IVIF_data[arr]=np.hstack((IVIF_data[arr],vals[arr])) #this assumes that IV is first
        return IVIF_data

    def summary_measure(self,celltrace,param):
        PSPsamples=np.zeros(len(self.sample_times))
        #PSPsamples[0]=param['meanpre']
        for index in range(len(self.sample_times)):
            norm_index=np.where(celltrace['psptime']/SEC_PER_MIN>self.sample_times[index])
            if len(norm_index[0]):
                psp_index=np.min(norm_index)
                PSPsamples[index]=np.nanmean(celltrace['normPSP'][psp_index-int(self.window/2):psp_index+int(self.window/2)])
            else:
                print ("calculating time samples: recording of ",param['exper'], "is shorter than ", self.sample_times[index], 'min')
                PSPsamples[index]=np.nan
        return PSPsamples
    ############################# FIXME: Edit this according to how parameters are saved in PatchAnal ##############################
    ##### Once induction parameters are saved, add them
    ##### possibly save Stim_interval differently
    def extract_params(self,exper_params,hs,data):
        cell_param={} #separate out parameter for just one headstage
        cell_param['exper']=exper_params['exper']+'_'+hs
        cell_param['date']='20'+exper_params['exper'].split('_')[0]
        for k,v in exper_params.items():
            if isinstance(v,dict):
                if hs in v.keys():                     
                    cell_param[k]=v[hs]
                elif k=='Stim_interval':
                    cell_param[k]=v
                elif k=='induction': 
                    if exper_params[k]['headstage']==hs: 
                        cell_param[k]=exper_params[k]['induct_inject'][0]
                    else:
                        cell_param[k]=0
                else:
                    print('unknown dictionary of experiment parameters', k,v)
            elif '(HS#' in k:  #pipette parameters
                hs_num=k.split('(HS#')[-1].split(')')[0]
                if hs_num in hs:
                    cell_param[k.split('(HS#')[0]]=v
            elif k in self.single_params:
                cell_param[k]=v
        for k,v in data.items():
            if isinstance(v,dict):
                if hs in v.keys():
                    cell_param[k]=v[hs]
        return cell_param

    def group_data(self,xvar,yvar): #FIXME: use this for IO curve, but xvar is stim current - so don't use SEC_PER_MIN
        self.grp_data = self.whole_df.groupby(self.sepvarlist)
        self.avg_PSP={}
        self.stderr_PSP={}
        self.samples={}
        self.minutes={}
        self.filenm={}
        for grp in self.grp_data.groups.keys(): #consider changing tuple to string using join
            self.filenm[grp]=self.construct_filename(grp)
            avg,std,count=grp_utl.exp_avg(self.grp_data.get_group(grp)[yvar])
            self.avg_PSP[grp]=avg
            self.stderr_PSP[grp]=std/np.sqrt(len(self.grp_data.get_group(grp)))
            self.samples[grp]=count
            self.minutes[grp]=grp_utl.exp_avg(self.grp_data.get_group(grp)[xvar])[0]/SEC_PER_MIN
    
    ########################### FIXME: Edit this according to exclusion criteria that are adopted ################
    def plot_bad(self):
        if len(self.BAD):
            fig,axes=pyplot.subplots(2,1)
            fig.canvas.manager.set_window_title('Problems')
            axes[0].set_ylabel('PSP baseline-bad slope or psp amp')
            axes[1].set_ylabel('PSP-nans')
            axes[1].set_xlabel('Time (min)')
            for exper,bad_data in self.BAD.items(): 
                time=bad_data['psptime']/SEC_PER_MIN
                if np.isnan(bad_data['normPSP']).any():
                    ### instead of plotting here, should put in separate container to plot later
                    print ("########## np.nan detected", exper)
                    p=axes[1].plot(time,bad_data['normPSP'],'+', label=exper)
                    color = p[-1].get_color()
                    nan_index=np.argwhere(np.isnan(bad_data['normPSP']))
                    axes[1].plot(time[nan_index],np.ones((len(nan_index))),'o', color=color)
                elif bad_data['slope']>self.slope_threshold or (np.abs(bad_data['slope'])-self.slope_std_factor*bad_data['slope_std']>0):
                    print("bad baseline", exper, "slope",round( bad_data['slope'],5), "+/-", round(bad_data['slopestd'],5))
                    axes[0].plot(time,bad_data['normPSP'],'.', label=exper)
                else:
                    print("PSP amp too low or insufficient traces", exper,'numnan=',sum(np.isnan(bad_data['normPSP'])))
                    axes[0].plot(time,bad_data['normPSP'],'x', label=exper)
                if bad_data['normPSP'][-1]==0.0:
                        print ("@@@@@@@@@ normPSP is 0, check PatchAnal for", exper)
                    #print 'OK: {}'.format(exper_param)
                print(bad_data['params'])       
            axes[0].legend(fontsize=8, loc='best')
            axes[1].legend(fontsize=8, loc='best')
            fig.canvas.draw()

    def write_stat_data(self):
        SASoutput = self.whole_df[self.single_params] #in to a 2d array you write it into SASoutput. 
        SASheader= '   '.join(self.single_params)
        SASoutput=np.column_stack((SASoutput,round(self.whole_df.meanpre,5)))
        for col in range(len(self.whole_df.PSPsamples[0])):
            pspmean=[round(row[col],5) for row in self.whole_df.PSPsamples]
            if not np.all([np.isnan(k) for k in pspmean]):
                SASoutput=np.column_stack((SASoutput,pspmean))
                SASheader += ' normPSP_'+str(self.sample_times[col])
        f=open("PARAMSforSAS.txt", 'w')
        f.write(SASheader +"\n")
        np.savetxt(f, SASoutput, fmt='%s', delimiter='   ')
        f.close()

    def group_to_word(self,grp):
        if isinstance(grp,tuple) or isinstance(grp,list):
            grp_nm='_'.join(grp)
        elif isinstance(grp,str):
            grp_nm=grp
        else:
            print('unknown group structure.  Not tuple or string or list')
            grp_nm=''
        return grp_nm

    def bar_graph_data(self,exclude_name): #this only plots and writes the 1st time sample
        lines=[]
        for grp in self.grp_data.groups.keys():
            grp_nm=[self.group_to_word(grp)]
            PSP_mean=self.grp_data.get_group(grp).PSPsamples.values
            line=grp_nm+[round(np.mean(PSP_mean,axis=0)[1]*100,5)]+[round((np.std(PSP_mean,axis=0)[1]/np.sqrt(len(PSP_mean)))*100,5)]
            lines.append(line)
            #print(line,np.std(PSP_mean,axis=0)[1])
        prefix=self.common_filename()
        f=open(prefix+"_BarGraphmeans.txt", 'w')
        f.write('  '.join(self.sepvarlist)+'    mean    sterr \n')
        np.savetxt(f, lines, fmt='%s', delimiter='   ')
        f.close()
        for grp in self.grp_data.groups.keys():
            PSP_mean=self.grp_data.get_group(grp).PSPsamples.values
            data_column=[ps[1]*100 for ps in PSP_mean] #convert the 1st time sample to percent
            f=open(self.filenm[grp]+"_points.txt",'w')
            if isinstance(grp,tuple):
                grp_vars=[g for i,g in enumerate(grp) if self.sepvarlist[i] not in exclude_name]
                columnheading=[prefix+'_'.join(grp_vars)]     #assumes you don't want theta as part of column name ; assumes drug comes before sex in ARGS           
            else:
                columnheading=[prefix+grp]
            f.write('_'.join(columnheading)+"_data\n")
            np.savetxt(f,data_column, fmt='%7.5f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
            f.close()
    
    def write_traces(self): #FIXME: use this for writing IO curve also?  Use current stim instead of minutes
        frac_2_percent = 100
        for gnum in self.grp_data.groups.keys():
            f=open(self.filenm[gnum]+".txt",'w')                      
            outputdata=np.column_stack((self.minutes[gnum], self.samples[gnum],
                                frac_2_percent*self.avg_PSP[gnum], frac_2_percent*self.stderr_PSP[gnum]))
            header="time "+self.filenm[gnum]+"count "+self.filenm[gnum]+"normPSP_AVG "+self.filenm[gnum]+"normPSP_SE "
            f.write(header+"\n")
            np.savetxt(f,outputdata, fmt='%7.5f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
            f.close()
    
    def write_IVIF(self):
        import scipy.stats
        self.IVIF_variables.remove('Im')
        for gnum in self.grp_data.groups.keys():
            f=open(self.filenm[gnum]+"_IVIF.txt",'w') 
            outputdata=self.grp_data.get_group(gnum).Im.mean()*1e12 #Convert to pA for printing, np.mean(self.grp_data.get_group(gnum).Im.values,axis=0)
            header='Im(pA)   '
            for yvar in self.IVIF_variables:
                #header=header+self.filenm[gnum]+yvar+'_avg   '+self.filenm[gnum]+yvar+'_sem  ' #longer name
                header=header+gnum+yvar+'_avg   '+self.filenm[gnum]+yvar+'_sem  '
                outputdata=np.column_stack((outputdata,np.nanmean(np.vstack(self.grp_data.get_group(gnum)[yvar].values),axis=0)))
                outputdata=np.column_stack((outputdata,scipy.stats.sem(self.grp_data.get_group(gnum)[yvar].values,axis=0,nan_policy='omit')))
            f.write(header+"\n")
            np.savetxt(f,outputdata, fmt='%7.5f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
            f.close()
            
    def common_filename(self):
        filnm=''
        common_params=[a for a in self.fname_vars]
        for p in common_params:
            if self.params.__getattribute__(p) is not None and p not in self.sepvarlist:
                filnm=filnm+str(self.params.__getattribute__(p))+'_'
        return filnm

    def construct_filename(self,paramgrp): 
        filnm=self.common_filename()
        self.common_filnm=filnm[0:-1]
        if len(self.sepvarlist)==1:
           paramgrp=[paramgrp]
        filnm=filnm+'_'.join([sv+str(v) for sv,v in zip(self.sepvarlist,paramgrp)])
        return filnm      

if __name__ =='__main__':        
    ARGS = "-sepvarlist region -plot_ctrl 111"      
    exclude_name=[] #['theta'] #use to exclude variable(s) from column name in _points files	        
    try:
        commandline = ARGS.split() #in python: define space-separated ARGS string
        do_exit = False
    except NameError: #if you are not in python, read in filename and other parameters
        commandline = sys.argv[1:]
        do_exit = True
	#### Age, drug, region, genotype allows you to select subset of data ###
    #### slope_std is value of exclusion criteria
    params=argp.ArgParserPatch(commandline,do_exit,1)
	#print('params=',params,'commandline =',commandline)            
    grp=GrpPatch(params)
    grp.pattern = grp.subdir+'*.npz'
    grp.outfnames = sorted(glob.glob(grp.pattern))
    if not len(grp.outfnames):
        sys.exit('************ No files found *********')
    if len(grp.outfnames):
        grp.read_data()
        grp.plot_bad()
        grp.group_data('psptime','normPSP') 
        if int(params.plot_ctrl[1])>0:
            plot_cols=int(params.plot_ctrl[1])
        else:
            plot_cols=None
        if int(params.plot_ctrl[0]):
            fig=grp_utl.plot_groups(grp.avg_PSP,grp.stderr_PSP,grp.minutes,grp.samples,grp.common_filnm,grp.sepvarlist,plot_cols)
        grp.write_traces() 
        grp.write_stat_data() 
        grp.bar_graph_data(exclude_name)
        grp.write_IVIF() 
        if int(params.plot_ctrl[2]):
            grp_utl.plot_corr(grp,['age'],'PSPsamples') ######### DEBUG
            grp_utl.plot_IVIF(grp,grp.IVIF_variables)
    #
    ########## NEXT STEPS: ################
    # Read in file to determine sex/pellet from ID 
    # Surgery date - calculate time since surgery
    #**** extract single rheobase and max_latency to use in corr plots?  should be only single value once valid experiments
    # extract Number of spikes during induction to use in corr plots?  After saving in patch Anal
    # Group IO data (after processing in patchAnal)
    ### IF needed, can add back in newcolumn_name - to take care of drug concentration - from GrpAvgPopSpikeClass
              
                
                
                
                
                
                
                
                