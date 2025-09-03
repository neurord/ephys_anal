# -*- coding: utf-8 -*-
"""
Created on Thu Jan 28 14:34:30 2021

@author: nehas
"""
import numpy as np
import sys  
import pandas as pd 
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 100)
import glob
import os
print (os.getcwd())
from matplotlib import pyplot
import ArgParser as argp
import PatchAnalPlots as pu
import GrpPlotUtil2 as grp_utl
import datetime
from PatchAnalPlots import line
from scipy import optimize

SEC_PER_MIN=60 #global
SYN_ROUTINE='Synaptic_StimDigCC'

class GrpPatch:
    def __init__(self, params): 
        self.params = params
        self.IDfile=params.IDfile
        self.subdir = params.outputdir
        self.slope_std_factor=params.slope_std 
        self.sepvarlist= params.sepvarlist
        self.sample_times=params.samp_time
        self.window=params.window
        self.plot_individuals=int(params.plot_ctrl[1])  #to plot PSP vs time for each experiment in a group
        self.plot_corr=int(params.plot_ctrl[2]) #to plot correlation between LTP and age or baseline epsp
        #additional parameters.  FIXME: add to arg parser
        self.minimum_sweeps=10#20 #5 min pre and 15 min follow-up
        self.slope_threshold=params.slope_thresh #fraction of change per sec.  Same as .0012 /minute or .036 in 30 min.
        self.nan_threshold=10 
        self.baseline_min=0#0.4
        self.baseline_max=2
        self.print_info=1 
        self.single_params=['region','genotype','age','drug','ID','time_to_induct','pre_num']
        self.fname_vars=['region','genotype','sex','drug']
        self.IVIF_variables=['Im','Vm','latency','num_spikes'] #also, APheight, APwidth, AHP_amp, AHP_time

    def dates(self,datestring,separator):
        if separator:
            date=datestring.split(separator)[0]
        else:
            date=datestring
        y=int('20'+date[0:2])
        m=int(date[2:4])
        d=int(date[4:6])
        newdate=datetime.date(y,m,d)
        return newdate

    def read_data(self): 
        self.BAD = {}
        DATAS = []
        PARAMS = []
        ANAL_PARAMS=[]
        IVIFset=[]
        IOset=[]
        for outfname in self.outfnames:
            with open(outfname, 'rb') as f:
                datadict = np.load(f,allow_pickle=True)
                data=datadict['data'].item()
                exper_param=datadict['params'].item()
                if exper_param['exper']=='250812_2':
                    exper_param['ID']='250530-B7-F2L' #not 250529-B9-F2L, wrong animal ID entered in metadata during experiment
                ## calculate time since surgery (females only), and age if necessary
                exper_date=self.dates(exper_param['exper'],'_')
                if exper_param['SxDate']: #no surgery date for males:
                    sx_date=self.dates(exper_param['SxDate'],'')
                    exper_param['sx_time']=(exper_date-sx_date).days
                else:
                    exper_param['sx_time']=np.nan #FIXME: might need to use something else to make correlation plots work
                if exper_param['age'] is None:
                    birth_date=self.dates(exper_param['ID'],'-')
                    exper_param['age']=(exper_date-birth_date).days                
                traces=datadict['trace'].item()
                IV_IF=datadict['IV_IF'].item()
                anal_params=datadict['anal_params'].item() #use if need to re-analyze experiment
                if 'IO' in datadict.keys():
                    IO=datadict['IO'].item()
                else:
                    IO={'amp':{'H2':[],'H1':[]}} #empty dictionary
                celltype=exper_param['celltype']
                print_params={k:v for k,v in exper_param.items() if not isinstance(v,dict) or k=='celltype'}
            if self.print_info and 'slope' in data.keys():
                print ("file read:", print_params,", baseline slope=", [round(sl,6) for sl in data['slope'].values()])
            ########## identify experiments that do not meet includsion criteria ##########
            for hs in celltype.keys():
                numnan=sum(np.isnan(traces['normPSP'][hs]))
                cell_param=self.extract_params(exper_param,hs,data)
                celltrace,induct_index,induct_index_avg=self.average_samples(traces,hs,exper_param['Stim_interval'])
                celltrace['PSPsamples']=self.summary_measure(celltrace,cell_param)
                #### Re-do slope to get intercept - not stored in .npz file. start slope at 10 min
                Aopt,Bopt,Bstd=self.new_slope(traces['amp'][hs],traces['psptime'][hs],10,induct_index)
                cell_param['Intercept']=Aopt
                #### Recalculate slope using entire baseline
                Aopt,Bopt,Bstd=self.new_slope(traces['amp'][hs],traces['psptime'][hs],0,induct_index)
                cell_param['slope10']=Bopt;cell_param['slope10_std']=Bstd;cell_param['Intercept10']=Aopt
                exclude = ((len(traces['normPSP'][hs])<self.minimum_sweeps) or
                                (np.abs(data['slope'][hs])>self.slope_threshold) or
                                (np.abs(data['slope'][hs])-self.slope_std_factor*data['slope_std'][hs]>0) or
                                #### Use slope determined from normPSP averaged over 1 minute
                                #(np.abs(cell_param['slope10'])>self.slope_threshold) or
                                #(np.abs(cell_param['slope10'])-self.slope_std_factor*cell_param['slope10_std']>0) or
                                numnan>self.nan_threshold or #do not use trace if too many missing popspikes
                                np.mean(data['meanpre'][hs])<self.baseline_min) #do not use if baseline psp is too small
                #If slope using all 10 minutes is bad, check the 5 min slope
                '''if (np.abs(cell_param['slope10'])>self.slope_threshold) or \
                        (np.abs(cell_param['slope10'])-self.slope_std_factor*cell_param['slope10_std']>0):
                    if (np.abs(data['slope'][hs])>self.slope_threshold) or \
                                (np.abs(data['slope'][hs])-self.slope_std_factor*data['slope_std'][hs]>0):
                        exclude=True'''
                if exclude: #print information message as to why experiment begin excluded
                    if len(traces['normPSP'][hs])<self.minimum_sweeps:
                        print ("!!!NOT ENOUGH goodtraces", exper_param['exper'], hs, len(traces['normPSP'][hs]))
                    elif data['meanpre'][hs]<self.baseline_min:
                        print ("!!!PSP amp too low", exper_param['exper'], hs, data['meanpre'][hs])
                    elif numnan>self.nan_threshold:
                        print('!!! Too many Nans',exper_param['exper'], hs, numnan)
                    elif np.abs(cell_param['slope'])>self.slope_threshold or (np.abs(cell_param['slope'])-self.slope_std_factor*cell_param['slope_std']>0):
                        print ("!!! BAD 5 min baseline", exper_param['exper'], hs, "5 min slope u,s", round(cell_param['slope'],7),round(cell_param['slope_std'],7))
                        if np.abs(cell_param['slope10'])>self.slope_threshold or (np.abs(cell_param['slope10'])-self.slope_std_factor*cell_param['slope10_std']>0):
                            print ("BAD 10 min baseline also: ", round(cell_param['slope10'],7),round(cell_param['slope10_std'],7))
                        else:
                            print(" 10 min baseline is OK:", round(cell_param['slope10'],7),round(cell_param['slope10_std'],7))
                    #add experiment to BAD to plot later
                    self.BAD[exper_param['exper']+'_'+hs]={'normPSP':celltrace['normPSP'],'psptime':celltrace['psptime'], 'meanpre':data['meanpre'][hs],
                                                        'params':{'celltype':celltype[hs],'age':exper_param['age'],'drug':exper_param['drug'],'region':exper_param['region']},
                                                        'slope':cell_param['slope'],'slope_std':cell_param['slope_std'], 'Intercept':cell_param['Intercept'],
                                                        'slope10':Bopt,'slope10_std':Bstd,'Intercept10':Aopt}
                    if self.print_info:
                        print ("***** excluding:", exper_param['exper']+'_'+hs)
                    next 
                else: #data meets all criteria
                    if np.isnan(traces['normPSP'][hs]).any(): #inform about nans, even if using data
                        print ("########## np.nan detected", numnan,'times in',exper_param['exper']+'_'+hs)
                    if traces['normPSP'][hs][-1]==0.0: #last response is 0, why????
                        print ("@@@@@@@@@ check PatchAnal for", exper_param['exper']+'_'+hs)
                    PARAMS.append(cell_param)
                    DATAS.append(celltrace)                    
                    IVIFset.append(self.IVIFdata(IV_IF,hs,cell_param))
                    IOset.append({'IOamp':IO['amp'][hs]}) 
                    ANAL_PARAMS.append(anal_params)       
        dfdata = pd.DataFrame(DATAS)
        dfparams=pd.DataFrame(PARAMS)
        df_ivif=pd.DataFrame(IVIFset)
        df_io=pd.DataFrame(IOset)
        df_anal=pd.DataFrame(ANAL_PARAMS)
        self.whole_df=pd.concat([dfparams,dfdata,df_ivif,df_io,df_anal], axis = 1)
        self.single_params.append('celltype')
        self.single_params.remove('pre_num')

    def ignore(self):
        ############ Select subset of files based on user specified criteria ##########  
        if self.params.sex:
            self.whole_df = self.whole_df[self.whole_df.Status == self.params.sex]
            print('ignoring experiments with sex not',self.params.sex)
        if self.params.age:
            self.whole_df = self.whole_df[self.whole_df.age >= self.params.age]
            print('ignoring experiments with age less than', self.params.age)
        if self.params.maxage:
            self.whole_df = self.whole_df[self.whole_df.age <= self.params.maxage]
            print('ignoring experiments with age greater than',self.params.age)
        if self.params.region:
            self.whole_df = self.whole_df[self.whole_df.region == self.params.region]
            print('ignoring experiments with region not ',self.params.region)
        if self.params.drug:
            self.whole_df = self.whole_df[self.whole_df.drug == self.params.drug]
            print('ignoring experiments with drug not ',self.params.drug)

    def new_slope(self,normPSP,time,starti,induct_index):
        pre_traces=normPSP[starti:induct_index]
        pre_time=time[starti:induct_index]
        no_nan_indices=~np.isnan(pre_traces)
        popt,pcov=optimize.curve_fit(line,pre_time[no_nan_indices],pre_traces[no_nan_indices])
        Aopt,Bopt=popt
        Astd,Bstd=np.sqrt(np.diag(pcov))  
        return Aopt,Bopt,Bstd
    
    def average_samples(self, traces, hs,stim_interval):
        time_samples={key:[] for key in traces.keys()}
        if isinstance(stim_interval,dict):
            follow_up=[v for r,v in stim_interval.items() if SYN_ROUTINE in r]
            interval=round(np.mean(follow_up))
        else:
            interval=stim_interval
        tracesPerMinute=int(SEC_PER_MIN/interval)
        calc_interval=np.diff(traces['psptime'][hs])
        induction_index=np.where(calc_interval>1.1*interval)[0][0]+1
        if induction_index%tracesPerMinute !=0:
            print ('**********baseline traces stopped without completing the minute - average will skip next samples - PROBLEM - fix code')
            #e.g. for samp in range(0,induction_index,tracesPerMinute)
            #THEN, for samp in range(induction_index,len(tr[hs]),tracesPerMinute)
        #DO NOT AVERAGE ACROSS INDUCTION INDEX
        for k, tr in traces.items():
            for samp in range(0,len(tr[hs]),tracesPerMinute):
                if samp<induction_index and samp+tracesPerMinute>induction_index:
                    time_samples[k].append(np.nanmean(tr[hs][samp:induction_index]))                    
                else:
                    time_samples[k].append(np.nanmean(tr[hs][samp:samp+tracesPerMinute]))
        return {k:np.array(v) for k,v in time_samples.items()}, induction_index, int(induction_index/tracesPerMinute)

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
                IVIF_data['max_latency']=np.nanmax(IVIF_data['latency'])
                if np.max(IVIF_data['num_spikes'])>0:
                    IVIF_data['rheobase']=[v for v in params['rheobase'].values()][0]
                else:
                    IVIF_data['rheobase']=np.nan
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
                if sum(np.isnan(bad_data['normPSP']))>self.nan_threshold:
                    ### instead of plotting here, should put in separate container to plot later
                    print ("########## np.nan detected", exper)
                    p=axes[1].plot(time,bad_data['normPSP'],'+', label=exper)
                    color = p[-1].get_color()
                    nan_index=np.argwhere(np.isnan(bad_data['normPSP']))
                    axes[1].plot(time[nan_index],np.ones((len(nan_index))),'o', color=color)
                elif np.abs(bad_data['slope'])>self.slope_threshold or (np.abs(bad_data['slope'])-self.slope_std_factor*bad_data['slope_std']>0):
                    print("bad baseline", exper, "slope",round( bad_data['slope'],7), "+/-", round(bad_data['slope_std'],7))
                    p=axes[0].plot(time,bad_data['normPSP'],'-', label=exper)
                    color = p[-1].get_color()
                    axes[0].plot(time,bad_data['Intercept']/bad_data['meanpre']+bad_data['slope']*(SEC_PER_MIN/bad_data['meanpre'])*time,linestyle=':', color=color)
                    axes[0].plot(time,bad_data['Intercept10']/bad_data['meanpre']+bad_data['slope10']*(SEC_PER_MIN/bad_data['meanpre'])*time,linestyle='-.', color=color)
                else:
                    print("PSP amp too low or insufficient traces", exper,'numnan=',sum(np.isnan(bad_data['normPSP']),'psp amp=', bad_data['meanpre']))
                    p=axes[0].plot(time,bad_data['normPSP'],'--', label=exper)
                if bad_data['normPSP'][-1]==0.0:
                    print ("@@@@@@@@@ normPSP is 0, check PatchAnal for", exper)
                    #print 'OK: {}'.format(exper_param)
                print(bad_data['params'])
            axes[0].set_title('dotted: 5 min slope, dot-dash: 10 min slope')       
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
                SASheader += 'baseline normPSP_'+str(self.sample_times[col])
        f=open(self.subdir+"PARAMSforSAS.txt", 'w')
        f.write(SASheader +"\n")
        np.savetxt(f, SASoutput, fmt='%s', delimiter='   ')
        f.close()

    def bar_graph_data(self,exclude_name): #this only plots and writes the 1st time sample
        lines=[]
        for grp in self.grp_data.groups.keys():
            grp_nm=[grp_utl.group_to_word(grp)]
            PSP_mean=self.grp_data.get_group(grp).PSPsamples.values
            line=grp_nm+[round(np.mean(PSP_mean,axis=0)[1]*100,5)]+[round((np.std(PSP_mean,axis=0)[1]/np.sqrt(len(PSP_mean)))*100,5)]
            lines.append(line)
            #print(line,np.std(PSP_mean,axis=0)[1])
        prefix=self.common_filename()
        f=open(self.subdir+prefix+"_BarGraphmeans.txt", 'w')
        f.write('  '.join(self.sepvarlist)+'    mean    sterr \n')
        np.savetxt(f, lines, fmt='%s', delimiter='   ')
        f.close()
        for grp in self.grp_data.groups.keys():
            PSP_mean=self.grp_data.get_group(grp).PSPsamples.values
            data_column=[ps[1]*100 for ps in PSP_mean] #convert the 1st time sample to percent
            f=open(self.subdir+self.filenm[grp]+"_points.txt",'w')
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
            f=open(self.subdir+self.filenm[gnum]+".txt",'w')                      
            outputdata=np.column_stack((self.minutes[gnum], self.samples[gnum],
                                frac_2_percent*self.avg_PSP[gnum], frac_2_percent*self.stderr_PSP[gnum]))
            header="time "+self.filenm[gnum]+"count "+self.filenm[gnum]+"normPSP_AVG "+self.filenm[gnum]+"normPSP_SE "
            f.write(header+"\n")
            np.savetxt(f,outputdata, fmt='%7.5f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
            f.close()
    
    def write_IVIF(self, ivif_dict):
        plot_vars=[a for a in self.IVIF_variables if a != 'Im']
        for group in ivif_dict.keys():
            f=open(self.subdir+self.filenm[group]+"_IVIF.txt",'w') 
            header='Im(pA)   '
            outputdata=ivif_dict[group]['Im']
            for yvar in plot_vars:
                header=header+self.filenm[group]+yvar+'_avg   '+self.filenm[group]+yvar+'_sem  '
                outputdata=np.column_stack((outputdata,ivif_dict[group][yvar],ivif_dict[group][yvar+'_ste']))
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
    ARGS = "Surgery_record -plot_ctrl 111"      #-sex FC -age 75
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
        grp_utl.read_IDfile(grp,'Sample',['Status']) #'Sample' has the unique identifier, ['Status'] is list of independent variables, e.g., sex, genotype to be added to df
        grp.ignore()
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
        ivif_dict=grp_utl.cluster_IVIF(grp,'Im',[a for a in grp.IVIF_variables if a != 'Im'],conversion=1e12 ) #convert to pA
        io_dict=grp_utl.cluster_IVIF(grp,'IOrange',['IOamp'],eps=.005)  #Use eps=.001-0.009 - since only specify two digits for stim) 
        grp.write_IVIF(ivif_dict) #FIXME: edit this to write IO curves
        if int(params.plot_ctrl[2]):
            grp_utl.plot_corr(grp,['age'],'PSPsamples') ######### FIXME: DEBUG
            grp_utl.plot_IVIF(grp,grp.IVIF_variables) #FIXME: edit this to write IO curves
            grp_utl.plot_IVIF_mean(grp,ivif_dict) #FIXME: edit this to write IO curves
    #
    ########## NEXT STEPS: ################
    # 2.Group IO data (after processing in patchAnal)
    #   
    # 3. extract Number of spikes during induction to use in corr plots?  After saving in patch Anal
    ### IF needed, can add back in newcolumn_name - to take care of drug concentration - from GrpAvgPopSpikeClass

                
                
                
                