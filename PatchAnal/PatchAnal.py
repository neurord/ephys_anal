

import h5py as h5
import numpy as np
import sys
import ArgParser as argp
from scipy.signal import find_peaks
import PatchAnalPlots as pu
import spike_utilities as su

HEADSTAGE_V={'H1':'_S2_', 'H2':'_S4_'}
HV_VAL_KEY={v:k for k,v in HEADSTAGE_V.items()}
HEADSTAGE_I={'H1':'_S1_', 'H2':'_S3_'}
CURRENT_THRESH=5e-12 #detect current injection with 5 picoAmp change in channels S1 or S3
SEC_PER_MSEC=0.001 #global
capacitive_artifact_points=2 
bad_psp=-25e-3 #above this value, probably an AP, not a PSP
MEASURES=['Timer_Time','Seg_start_time'] #code assumes that Seg_start_time indicates time of digital stimulation
stim_per_burst=4 #FIXME: for theta burst - different value for LTD - add to ArgParser

class PatchAnal():
    def __init__(self,params): #change to PopSpike(args)
        self.datadir = params.datadir
        self.outputdir = params.outputdir
        self.filename_ending=params.filetype
        self.experiment=params.experiment
        self.region=params.region
        self.ID=params.ID
        self.genotype=params.genotype
        self.age=params.age
        self.drug=params.drug
        self.graphs=params.graphs
        self.plotstart=params.plotstart 
        self.slope_std_factor=params.slope_std
        self.ss_dur=params.ss_dur
        self.PSPstart=params.PSPstart+params.decay #Extracted from notebook file if one exists
        self.basestart=params.basestart
        self.base_dur=params.base_dur
        self.window=params.window #filter window for finding peak PSP
        self.APthresh=params.APthresh  
        self.thresheight=params.thresheight
        self.refract=params.refract
        self.max_risetime=params.max_risetime
        self.min_risetime=params.min_risetime
        self.threshval=params.threshval
        self.headstages=params.headstages
        self.IOrange=np.array(params.IOrange)
        self.digstim=params.digstim 
        self.Stim_interval=params.PSP_interval
        self.artifact_decay=params.decay #time for decay of stimulation artifact
        self.induction=params.induction
        self.baseline_time=params.base_time
        if len(params.headstages)==len(params.celltype): 
            self.celltype={h:params.celltype[i] for i,h in enumerate(self.headstages)}
        else:
            print('*********** Need to specify celltype for each valid headstage ***********')
            exit()
        self.params={'exper':self.experiment,'region':self.region, 'genotype':self.genotype, 'age':self.age, 'drug':self.drug, 'ID':self.ID, 'celltype': self.celltype}
        self.anal_params={
                     'PSPstart':self.PSPstart, 'basestart': self.basestart, 'base_dur':self.base_dur, 'ss_dur': self.ss_dur, 'base_time': self.baseline_time,
                     'window': self.window, 'IOrange':np.array(self.IOrange),'digstim':self.digstim, 'threshval': self.threshval, 'decay':self.artifact_decay,
                     'APthresh':self.APthresh,'refract':self.refract,'max_risetime':self.max_risetime,'min_risetime':self.min_risetime}

    def read_datafile(self,question=True): #NOTE: create separate read_datafile if reading exported ibw files
        self.filename=self.datadir+self.experiment+self.filename_ending 
        self.data = h5.File(self.filename,"r") 
        #print(self.data.keys()) #['Data', 'Meta', 'Analysis', 'Routines', 'Paradigms', 'Images', 'Logging', 'Notebook', 'DataWaveNotes', 'ExperimentStructure']
        #self.data['Data'].keys() - list of all routines, e.g. R1_S1_IV_CC, which has shape (num_points,num_sweeps)
        #dimensions of data arrays is number of time points x number of traces
        #['Paradigms'] - keys are executed paradigms, each of which is text array giving the routines and paradigms
        #['Routines'] - keys are executed routines, each of which encoded
        #['Logging'] - list of paradigms and routines, start and stop times
        #['DataWaveNotes'] - #first row of each entry contains: units, number of points, min and max values, dt for traces

        self.routines=list(self.data['Data'].keys())
        psp_list=[r for r in self.routines if r.endswith('Synaptic_StimDigCC')] #FIXME: hardcoded name of baseline and followup
        #Determine which routines are pre theta and which are post-theta
        self.induction_list=[r for r in self.routines if self.induction in r] 
        if len(self.induction_list):
            first_burst=self.induction_list[0]
            self.induct_headstage='H'+str(int(first_burst[-1])) #identify which headstage had depol during theta
            self.params['pre_num']=int(first_burst.split('_')[0][1:])-1 #number of routine that has baseline PSPs is one before first theta burst
        elif len(psp_list):
            StimDigCC=[a.split('_')[0] for a in psp_list]
            Stim_routines=np.unique(StimDigCC)
            self.params['pre_num']=int(Stim_routines[0][1:])
        else:
            print('no psp or induction routines found')
        #
        self.psp_dict={}
        self.IV_IF_dict={}
        routine_list=[]
        for row in self.data['DataWaveNotes'].__array__(): #units, number of points, min and max values, dt for traces
            routine_list.append(row[0].decode('UTF8'))
        for p in psp_list: #dict of traces and dt info for psps.  Possibly create IV and IF dictionaries
            self.psp_dict[p]=routine_list.index(p)
            if p.startswith('R'+str(self.params['pre_num'])):
                self.pre_trace=p
        #
        IV_list=[r for r in self.routines if 'IV_CC' in r]#FIXME: hardcoded name of IV and IF curves
        IF_list=[r for r in self.routines if r.endswith('IF_CC')]#FIXME: hardcoded name of IV and IF curves
        #dict of traces and dt info for IV and IF 
        self.IV_IF_dict={p: routine_list.index(p) for p in IV_list+IF_list} 
        #
        self.induct_dict={p:routine_list.index(p) for p in self.induction_list}
        ##dict of traces and dt info for IO curves.
        IOtest=[r for r in self.routines if 'StimIOtest' in r]  #FIXME: hardcoded name of IO curves
        self.IO_dict={p:routine_list.index(p) for p in IOtest}
        #
        self.dt_dict=self.psp_dict|self.IV_IF_dict|self.induct_dict|self.IO_dict #used to convert from I_dt to V_dt
        print('PSP:', self.psp_dict,'\nTheta:', self.induct_dict)
        self.log_rows= self.data['Logging']['Logging'].__array__()[0].decode('UTF8').split('\r')
        #Questions: how to determine time of digital stimulation
        #initialize sweep_time, to be filled if notebook file exists
        routines=np.unique([r.split('__')[-1] for r in self.data['Routines'].keys()])
        self.sweep_time={r:[] for r in routines} #obtain from notebook file if available
        self.params['Stim_interval']={r:self.Stim_interval for r in routines} #initialize here, in case no notebook file

    def read_notebook(self): #extract routine and sweep time from notebook file, if it was saved
        #
        def parse_measures(line,measure): #measure is the name of the measurement we are looking for, e.g. Timer_time or Seg_start_time
            parts=line.split(',')
            measure_names={i:p for i,p in enumerate(parts) if measure in p} #could just be list?
            value=[float(parts[i+1].split()[0]) for i in measure_names.keys()]
            return value #empty list if line="" or measure no found

        import glob
        self.routine_time={} #initialized regardless of whether notebook file exists
        PSPstart_dict={} #if notebook file exists, ignore the PSPstart read in ArgParser
        #nfname=self.datadir+'Notebook_20'+self.experiment.split('_')[0]+'*.txt'  #find all Notebook files with correct date
        nfname=self.datadir+'Notebook_20'+'*.txt'  #find all Notebook files with correct date
        files=glob.glob(nfname)
        self.theta_in_notebook=False
        if len(files):
            findex,file_found=su.find_notebook_file(files,self.experiment)

        else:
            file_found=False
        if file_found: 
            print('Using notebook file:',files[findex])  
            with open(files[findex],'r') as myfile: 
                all_lines=myfile.readlines()
            self.series_start={}
            for i,line in enumerate(all_lines):
                if line.startswith('Series'): #series name
                    series=line.split(': ')[-1][0:-1]
                    self.series_start[series]=[i+1] #starting line number for series/routine name
                    if self.induction in series: 
                        self.theta_in_notebook=True
                elif line=='\n' and len(self.series_start[series])<2: #end of series/routine, except for last series. Requires that these lines are encountered AFTER Series
                    self.series_start[series].append(i) 
                elif line.startswith('1 - '): #remaining lines: read meta data
                    if 'Animal Age' in line:
                        self.age=self.params['age']=line.split(', ')[-1][0:-1]
                    elif 'Animal Identifier' in line:
                        self.ID=self.params['ID']=line.split(', ')[-1][0:-1]
                    elif 'Tissue Preparation' in line:
                        self.region=self.params['region']=line.split(', ')[-1][0:-1] #override default
                    elif 'Animal Preparation' in line:
                        self.SurgeryDate=self.params['SxDate']=line.split(', ')[-1][0:-1] #override default
                    elif 'Bath Solution' in line:
                        self.drug=self.params['drug']=line.split(', ')[-1][0:-1] #override default
                #else:
                    #do nothing with remaining lines
            if len(self.series_start[series])<2: #if last series doesn't have an end, add the last line
                self.series_start[series].append(i+1)
            if not self.ID:
                print('!!!!!!!!!!! Animal ID is not in Notebook file, please re-do analysis and enter ID, age and OVX date on command line !!!!')
                #possible system exit at this point
            #loop over the series/routines, and read Timer_time for all sweeps
            for series,[start,end] in self.series_start.items():
                for ln in range(start,end):
                    if all_lines[ln].startswith(series.split('_')[0]) and all_lines[ln+1].startswith(series.split('_')[0]): #make sure the line isn't broken in the middle
                        full_line=all_lines[ln]
                    elif all_lines[ln].startswith(series.split('_')[0]) and not all_lines[ln+1].startswith(series.split('_')[0]): #next line is continuation
                        full_line=all_lines[ln]+all_lines[ln+1]
                        ln+=1 #skip the next line
                    else:  #Possibly end is incorrect value (for last Series) or other incorrect lines have been inserted
                        full_line=""
                    if len(full_line):
                        t=parse_measures(full_line,MEASURES[0])
                        if len(t):
                            self.sweep_time[series].append(t[0]) #Timer_time
                        PSPstart=parse_measures(full_line,MEASURES[1]) #find Segment start time - when digital stimulation applied - may be 0, 1 or two values
                        if len(PSPstart)==0:
                            PSPstart_dict[series]=self.PSPstart  #if not found (PSPstart=[]) then use default value 
                        elif len(PSPstart)==1:
                            PSPstart_dict[series]=PSPstart[0] + self.artifact_decay #save single value
                        else:  
                            PSPstart_dict[series]=[x+self.artifact_decay for x in PSPstart] #list

                self.routine_time[series]=self.sweep_time[series][0] #time of routine start is ~Timer_Time of 1st sweep
                if len(self.sweep_time[series])>=2:
                    self.params['Stim_interval'][series]=np.mean(np.diff(self.sweep_time[series]))
                if series.startswith('R1_'):
                    self.breakin_time=self.sweep_time[series][0] #approximate breakin time as Timer_Time of 1st sweep of 1st routine - actual breakin is smaller than this
                    #### if timer reset is in log, then can define self.breakin_time as 0
            self.PSPstart=PSPstart_dict #once PSPstart for all routines is found, can override the default/argParser value
        else:
            print('Notebook file not found using pattern:',nfname)

    def time(self): #### Extract routine time from logging.  Needed if notebook file not found or not all routines in notebook file
        def hours_to_sec(tim): #convert from AM/PM to seconds
            parts=tim.split(' ')
            actual_time=parts[-2]
            h,m,s=actual_time.split(':')
            #NOTE, this does deal properly with midnight
            if parts[-1]=='PM' and h!='12':
                h=int(h)+12                        
            return sum(x * int(t) for x, t in zip([3600, 60, 1], [h,m,s])) 
        #obtain routine start from Logging
        start_rows=[]; exp_start=[];timer_reset=None
        for lr in self.log_rows:
            if "Routine acquisition started" in lr:
                start_rows.append(lr)
            if "Starting new Paradigm" in lr:
                exp_start.append(lr)
            if "Tag_Timer_Reset" in lr:
                timer_reset=lr
        if timer_reset:
            start_time=hours_to_sec(timer_reset.split('\t')[0].split(',')[-1]) #breakin_time
            self.breakin_time=0
        else:
            start_time=hours_to_sec(exp_start[0].split('\t')[0].split(',')[-1]) #approximate breakin_time as start time of 1st paradigm - breakin is probably earlier than this
            if not hasattr(self,'breakin_time'):
                self.breakin_time=start_time
            print('no timer reset, using Paradigm Start Time as breakin time for routine time from logging, and 1st sweep time as breakin time for routine time from notebook')
        print('start_time:', exp_start[0].split('\t')[0].split(',')[-1], ', in seconds',start_time, ', breakin time:',self.breakin_time)
        for sr in start_rows: #only use routine start from logging if routine time not available from notebook file
            key=sr.split('\t')[3]
            time=hours_to_sec(sr.split('\t')[0].split(',')[-1])
            if key in self.routine_time.keys():
                print('routine time for',key,' from log =',time-start_time,'vs Notebook =',self.routine_time[key]-self.breakin_time)
            else:
                self.routine_time[key]=time-start_time #comparable to Timer_Time? 
        #########################################################################
        #           If notebook file 
        #               use timer_time minus 1st sweep time as routine time
        #               breakin time = timer_time of 1st sweep  
        #               If timer-reset
        #                   breakin_time = 0
        #           Else (not notbooke file, or not all routines in notebook file):
        #               use routine time from logging, but subtract 1st paradigm time
        #               If timer_reset  
        #                   subtract timer_reset time
        #           If timer_time of induction is available (in notebook file)
        #               time_to_induct = timer_time of 1st induction trace - breakin_time, which is either 0 (if timer_reset) or time of 1st sweep (no reset)
        #           else (use routine time of induction from log)
        #                time_to_induct = routine_time - start_time
        #                do not subtract breakin_time (if timer_reset - value is 0, if no timer_reset - start_time has already been subtracted)
        #########################################################################

    def meta_data(self): #expand this once we can get animal ID, etc?
        pipette=[]
        for lr in self.log_rows:  #obtain pipette resistance, etc from Logging
            if '(HS#' in lr:
                pipette.append(lr)
        for pr in pipette: #add pipette resistance, etc to parameters
            info=pr.split('\t')[2].split(' = ')
            self.params[info[0]]=info[1]

    def init_arrays(self): #NOTE: separate out self.Vm_psp, and store traces, if need to read data from exported ibw
        self.Vm_psp={}; self.Vm_IV_IF={}; self.IO_psp={}
        self.num_traces={}
        self.IV_IF_spikes={}
        self.induction_time=[]
        for headstage in self.headstages: ### initialize dictionary of arrays - one dict entry for each headstage
            self.Vm_psp[headstage]= {r:i for r,i in self.psp_dict.items() if HEADSTAGE_V[headstage] in r} #list of routines for psp extract, and number with dt
            self.IO_psp[headstage]= {r:i for r,i in self.IO_dict.items() if HEADSTAGE_V[headstage] in r}
            for r in self.Vm_psp[headstage].keys():
                print(r,np.shape(self.data['Data'][r].__array__()))
            self.num_traces[headstage]=sum([np.shape(self.data['Data'][r].__array__())[-1] for r in self.Vm_psp[headstage].keys()])
            self.Vm_IV_IF[headstage]= {r:i for r,i in self.IV_IF_dict.items() if HEADSTAGE_V[headstage] in r}#list of routines for IV-IF analysis and number with dt
            self.IV_IF_spikes[headstage]={r:{} for r in self.IV_IF_dict.keys() if HEADSTAGE_V[headstage] in r} #will create separate recarrays for spikes in each burst
        self.IV_IF={h:{} for h in self.headstages} #will create recarrays for sweep level values
        self.max_latency={h:{} for h in self.headstages} #max time between current and AP across all sweeps
        self.rheobase={h:{} for h in self.headstages} #minimum current that produces AP
        #initialize arrays to hold psp values - possibly convert to single recarray
        self.psp={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages } #peak value of psp
        self.pspamp={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages } #psp amplitude - subract RMP from psp
        self.psptime={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages } #time that trace starts from start of baseline
        self.peaktime={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages } #time within trace of psp peak
        self.RMP={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages} #resting membrane potential
        self.hyper={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages} #hyperpolarization value (testing access resistance)
        self.Raccess={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages} #hyper/Iaccess
        self.dV_2ms={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages} #hyper/Iaccess
        self.normpsp={headstage:np.zeros(self.num_traces[headstage]) for headstage in self.headstages } #PSP normalized to baseline mean

        ####################################### Next Steps ###############################
        #2. time of dig stimulation within trace - segment time?? for theta and StimDigCC?
        #3. bad traces: for each trace, calculate access resistance - if too many are bad - exclude the experiment, here or in GrpAvg?

    def find_inject(self,r): #use this to find injection time during thetaburst, IV,IF
        #name of routine for current determined from name for voltage
        #timing must be same for all inject in this function
        I_r=self.name_inject(r)
        inject=self.data['Data'][I_r].__array__()#[:,0] 
        deriv=np.diff(inject,axis=0)
        current_changes0=find_peaks(np.abs(deriv[:,0]),CURRENT_THRESH) #find time of current inject from 1st trace
        timepoints=current_changes0[0]
        if np.shape(deriv)[-1]>1:
            current_changes1=find_peaks(np.abs(deriv[:,1]),CURRENT_THRESH) #find time of current inject from 2nd trace
            if len(current_changes1[0])>len(current_changes0[0]):
                print('**** PROBLEM: time of current injection is different for each trace for', r, '!!!', len(current_changes1[0]), 'vs', len(current_changes0[0]))
                timepoints=current_changes1[0]
        return timepoints,inject  #time points is the time just before the change in current, due to shift in derivative

    def name_inject(self,r):
        parts=r.split('_')
        channel='_'+parts[1]+'_'
        if channel in HV_VAL_KEY.keys():
            headstage=HV_VAL_KEY[channel]
        parts[1]=HEADSTAGE_I[headstage][1:-1]
        new_r='_'.join(parts)
        return new_r

    def get_dt(self,routine_num): #use this to find dt for any routine 
        #R1,R2, etc = number of routine
        #S1: channel 1 current,S2: channel 1 voltage;S3-channel2 current;S4-channel 2 current
        text=self.data['DataWaveNotes'].__array__()[routine_num][2].decode('UTF8')
        dt=float(text.split('XINC=')[-1].split(';')[0]) #needs to be identifical for both headstages and all psp routines
        return dt
    
    def Itime_to_Vtime(self,Itimepoints,Vdt,r):
        Idt=self.get_dt(self.dt_dict[self.name_inject(r)])
        factor=int(Idt/Vdt)
        if np.abs(Idt/Vdt-factor)>1e-12: #FIXME: replace with SMALL_NUMBER
            print('PROBLEM: Vdt and Idt are not related by integer factors!!!')
        Vtimepoints=[tp*factor for tp in Itimepoints] #Vtimepoints could be multiple points prior to current onset
        return Vtimepoints, factor

    def get_routine_start(self, r):
        parts=r.split('_')
        del parts[1] #parts[0] is routine number, parts[1] is channel, remainder is name
        routine_name='_'.join(parts)
        routine_start=self.routine_time[routine_name] #needed for summary plot and slope
        return routine_start,routine_name
    
    def analyze_theta(self):  #extract time of theta - to calculate time to induction
        def induct_params(routine_name):
            #determine burst freq, induct_freq and train_interval 
            if isinstance(self.PSPstart,dict):
                if isinstance(self.PSPstart[routine_name],list):
                    burst_interval= np.mean(np.diff(self.PSPstart[routine_name])) #20 ms for Theta, nan for 20 Hz
                    self.params['induction']['burst_interval']=burst_interval
                else:
                    self.params['induction']['burst_interval']=np.nan
            else:
                self.params['induction']['burst_interval']=np.nan #unknown
            self.params['induction']['induct_interval']=self.params['Stim_interval'][routine_name] # 95 ms or 1/20 Hz for LTD
            #The following is kluge.  Don't want to specify routine name (e.g. ThetaBurst), cuz don't know name of 20Hz or other possible induction protocols
            #So, instead exclude specific, non-induction protocols.  Will be problematic if someone adds a protocol
            #Could do something similar to MEASURES and define induction protocols at top
            routine_times=[x for nam,x in self.routine_time.items() if self.induction in nam]
            if len(routine_times):
                self.train_interval=np.mean(np.diff(routine_times))
                self.params['induction']['train_interval']=self.train_interval
            else:
                self.params['induction']['train_interval']=np.nan
                print('******* Unable to determine train interval - fix "routine_times=" in "analyze_theta" !!!!!!!!!!!!!!!!')
            return #not needed

        def find_stim_times(routine_name):
            adjust=-self.artifact_decay#self.min_risetime-self.artifact_decay #artifact decay is too long for labeling stim time, to show it occurs prior to AP
            if isinstance(self.PSPstart,dict):
                if isinstance(self.PSPstart[routine_name],list):  #Might need to use isinstance instead of length
                    #FIXME: APs seem to occur immediately after stim, or sometimes ~4 ms later
                    stim_time=np.array([self.PSPstart[routine_name][0]+adjust+i*np.mean(np.diff(self.PSPstart[routine_name])) for i in range(stim_per_burst)])
                else:
                    stim_time=np.array(self.PSPstart[routine_name])
            else:
                stim_time=[0.015,0.035,0.055,0.075] #FIXME: make this parameter in ArgParser, or make this empty list?
            return stim_time
        
        def ISI_anal(AP_times, stim_time):
            #if there is no burst frequency, just use time of single stim
            for pt in AP_times: #for each spike)
                delta=pt-stim_time #positive delta means spiketime AFTER psp - forward pairing
                minloc=np.abs(delta).argmin()
                self.min_interval.append(delta[minloc])
                if pt>stim_time[0] and pt<stim_time[-1]:
                    if self.min_interval[-1]>0:
                        self.other_interval.append(delta[minloc+1]) #AP closest to previous stim. oher interval time from AP to next stim
                    else:
                        self.other_interval.append(delta[minloc-1]) #AP closest to next stim. other interval time from AP to previous stim
            return

        self.num_spikes={r:{} for r in self.induct_dict.keys() if HEADSTAGE_V[self.induct_headstage] in r} #will create separate recarrays for spikes in each burst
        self.spikes={r:{} for r in self.num_spikes.keys()}
        self.params['induction']={'headstage':self.induct_headstage}
        self.min_interval=[]; self.other_interval=[] #store PSP to spike intervals
        for r in self.spikes.keys():
            timepoints,inject=self.find_inject(r) 
            self.induct_dt=self.get_dt(self.induct_dict[r])
            self.induction_time.append(self.get_routine_start(r)[0]) #use to calculate time from break-in to induction
            routine_start,routine_name=self.get_routine_start(r)
            stim_time=find_stim_times(routine_name) #should be the same for each routine
            if len(timepoints): #determine value of current injection during theta
                Vtimepoints,factor=self.Itime_to_Vtime(timepoints,self.induct_dt,r)
                self.depol_startpt=Vtimepoints[0]+factor #make sure AFTER current onset 
                self.depol_endpt=Vtimepoints[1]-capacitive_artifact_points #make sure capacitive artifact excluded - time BEFORE current offset
                self.params['induction']['induct_inject']=inject[timepoints[1]]-inject[0] # timpoint[1] is end of current injection
            else:
                print('no current injection during plasticity induction', r)
                self.params['induction']['induct_inject']=0
                self.depol_startpt=0 #replace this with some time relative to digital stimulation
                self.depol_endpt=len(self.data['Data'][r].__array__())
            #find spikes in each sweep
            for num in range(np.shape(self.data['Data'][r].__array__())[-1]):
                trace=self.data['Data'][r].__array__()[:,num]
                #Refractary period after stim artifact is preventing finding some APs
                peaks=find_peaks(trace[self.depol_startpt:self.depol_endpt],self.APthresh,width=10) #width in points, stim artifact interferes with spike detection if distance is used
                if len(peaks[0]): #verify AP and not noise, extract spike time and risetime
                    num_peaks,_,peak_times=self.find_spikes(self.induct_headstage,r,num,trace[self.depol_startpt:self.depol_endpt],peaks[0],self.induct_dt,routine_type='induct')
                else:
                    num_peaks=0
                    peak_times=[]
                self.num_spikes[r][num]=num_peaks #num spikes within burst
                if len(peak_times) and len(stim_time): 
                    ISI_anal(peak_times,stim_time) 

        induct_params(routine_name)
        if self.theta_in_notebook:
            print('timer time available for induction')
            #use timer time - breakin time for time to induction
            self.params['time_to_induct']=self.induction_time[0]-self.breakin_time #induction_time uses timer_time
        else:
            print('timer time NOT available for induction, using time of induction - time of 1st paradigm')
            #subtract 
            self.params['time_to_induct']=self.induction_time[0] #induction_time uses clock time converted to seconds, and then time of 1st paradigm is subtracted already
        return stim_time

    #once time of dig-stim is known, calculate time between AP and dig-stim?
    #possibly calculate AP frequency within burst?
    def find_spikes(self,h,r,s,wave,peaks,dt,routine_type): #FIXME: change to induct vs IV_IF?
        time=np.arange(len(wave))*dt
        riset=[];Vthresh=[];APtime=[];APpeak=[];APheight=[]
        delete_list=[]; sweep_freq=None
        for i,peakpt in enumerate(peaks):
            start=max(0,peakpt-int(self.max_risetime/dt))
            if start<len(wave)-2:
                y=wave[start:peakpt+2]
                x=time[start:peakpt+2]
                yderiv=np.diff(y)
                ythresh=self.threshval*yderiv.max()
                threshV=y[1:][yderiv>ythresh].min()
                threshtime=x[1:][yderiv>ythresh].min()
                risetime=time[peakpt]-threshtime
                if risetime<self.min_risetime or wave[peakpt]-threshV<self.thresheight: #(e.g. 10 or 50 mV)
                    delete_list.append(i)
                    #remove peakpt from peaks
                else:
                    APtime.append(time[peakpt]) #APtime is relative to start of inject
                    riset.append(risetime) #start time of AP is APtime-risetime
                    Vthresh.append(threshV)
                    APpeak.append(wave[peakpt])
                    APheight.append(wave[peakpt]-threshV)
        if len(delete_list):
            peaks=np.delete(peaks,delete_list)
        num_peaks=len(peaks)
        #if num_peaks==0 (all peaks are noise) don't store these
        if routine_type=='IVIF' and num_peaks:
            sweep_freq=np.mean(1/np.diff(APtime))
            width=su.spike_width(wave,time,peaks,APheight,Vthresh)
            ahp,ahp_pt=su.spike_ahp(wave,peaks,(np.array(width)/dt).astype(int),Vthresh)
            if len(ahp)<num_peaks: #redimension all arrays of AP characteristics to drop the last spike
                num_peaks=len(ahp)
                for AParray in [APtime,riset,Vthresh,APpeak,APheight,width]:
                    AParray.pop()
            ahp_time=time[ahp_pt]-APtime
            self.IV_IF_spikes[h][r][s]=np.rec.fromarrays((APtime,riset,Vthresh,APpeak,APheight,width,ahp,ahp_time),names='APtime,risetime,Vthresh,APpeak,APheight,APwidth,AHP_amp,AHP_t')
        elif routine_type=='induct' and num_peaks: #only store peak value, peak time and risetime for theta protocol
            self.spikes[r][s]=np.rec.fromarrays((APtime+self.depol_startpt*self.induct_dt,riset,APpeak),names='APtime,risetime,APpeak')
        elif num_peaks:
            print('!!!!!!!!!!!!!!!!! un-recognized routine type !!!!!!!!!!!!!!!!!!!!! ')
        else:
            print('no peaks found in', r,s)
        return num_peaks, sweep_freq, APtime

    def IV_IF_anal(self):
        def inject_time(timepoints,dt,r): #find injection current for both IV and IF
            #two current injections during PSP.  Use the last one for access resistance
            Vtimepoints,factor=self.Itime_to_Vtime(timepoints,self.IV_dt[r],r)
            inject_endpt=Vtimepoints[1]-capacitive_artifact_points #make sure capacitive artifact excluded
            inject_startpt=Vtimepoints[0]+factor  #make sure AFTER current onset. not used? 
            ss_startpt=inject_endpt-int(self.ss_dur/dt) #duration for measuring steadystate Vm
            inject_start=inject_startpt*dt
            base_end=Vtimepoints[0]-capacitive_artifact_points #is this early enough?  
            base_start=max(0,Vtimepoints[0]-int(self.basestart/dt)) 
            print('onset', Vtimepoints[0], round(inject_start,5),'s, offset',round(inject_endpt*dt,5))
            return inject_startpt,inject_endpt,base_start,base_end,ss_startpt
        
        def calc_mean(spikes_sweep,apdict,measure):
            new=[np.nan for i in range(len(spikes_sweep))]
            for i,sweep in apdict.items():
                new[i]=np.mean(sweep[measure])
            return new
        
        self.IV_dt={}
        onset={h:{} for h in self.headstages} #time within sweep that current starts and stops
        offset={h:{} for h in self.headstages}
        for headstage in self.Vm_IV_IF.keys():
            for r in self.Vm_IV_IF[headstage].keys():
                spikes_sweep=[];APfreq=[]
                self.IV_dt[r]=self.get_dt(self.IV_IF_dict[r])
                #above provides dt of Vm
                #determine time and value of hyperpolarizing pulses
                timepoints,inject=self.find_inject(r) #time when current starts and stops 
                inj_startpt,inj_endpt,base_start,base_end,ss_start=inject_time(timepoints,self.IV_dt[r],r) #time to measure ss before and during injection
                traces=self.data['Data'][r].__array__()
                Im=np.mean(inject[timepoints[0]:timepoints[1]],axis=0)-inject[0]  #timepoints[0] is just before inject, so cannot use that without adding small integer
                Vm=np.mean(traces[ss_start:inj_endpt],axis=0)-np.mean(traces[base_start:base_end],axis=0) # value of ss deltaV
                rect=su.rectification(traces,Im, inj_startpt,inj_endpt,self.window)-np.mean(traces[base_start:base_end],axis=0) #for negative current inject, determine whether rectification
                for num in range(np.shape(traces)[-1]):  #num indexes traces / sweeps
                    peaks=find_peaks(traces[inj_startpt:inj_endpt,num],self.APthresh,distance=int(self.refract/self.IV_dt[r])) #AP are defined as crossing zero, exceeding APthresh
                    if len(peaks[0]):
                        num_spikes,sweep_freq,_=self.find_spikes(headstage,r,num,traces[inj_startpt:inj_endpt,num],peaks[0],self.IV_dt[r],routine_type='IVIF')                        
                    else:
                        num_spikes=0
                        sweep_freq=np.nan
                    spikes_sweep.append(num_spikes)
                    APfreq.append(sweep_freq)
                #mean spike values per sweep
                mean_height=calc_mean(spikes_sweep,self.IV_IF_spikes[headstage][r],'APheight')
                mean_width=calc_mean(spikes_sweep,self.IV_IF_spikes[headstage][r],'APwidth')
                mean_ahp=calc_mean(spikes_sweep,self.IV_IF_spikes[headstage][r],'AHP_amp')
                mean_ahpt=calc_mean(spikes_sweep,self.IV_IF_spikes[headstage][r],'AHP_t')
                latency=[np.nan for i in range(len(spikes_sweep))]
                for i,sweep in self.IV_IF_spikes[headstage][r].items(): #use -1 instead of nan??
                    latency[i]=sweep.APtime[0]
                self.IV_IF[headstage][r]=np.rec.fromarrays((Im,Vm,rect,spikes_sweep,latency,mean_height,mean_width,mean_ahp,mean_ahpt,APfreq),
                                                           names='Im,Vm,rect,num_spikes,latency,APheight,APwidth,AHP_amp,AHP_time,APfreq')
                #values that are once per routine.  add in: IR per sweep=Vm/Im, frequency?  or calculate later
                if np.max(spikes_sweep)>0:
                    self.max_latency[headstage][r]=np.nanmax(latency)  #or if latency=-1, can just use np.max
                    lowest_spike=np.min(np.where(np.array(spikes_sweep)>0))
                    self.rheobase[headstage][r]=Im[lowest_spike]
                onset[headstage][r]=inj_startpt*self.IV_dt[r] 
                offset[headstage][r]=inj_endpt*self.IV_dt[r]
        self.params['IV_onset']=onset
        self.params['IV_offset']=offset

    def calc_psp(self):
            
        def access_monitor(timepoints,inject, dt,r): #modify this to find injection current for IV and IF
            #two current injections during PSP.  Use the last one for access resistance
            #above could be used to find time of inject during IV-IF curves - make independent function
            Vtimepoints,factor=self.Itime_to_Vtime(timepoints,dt,r)
            self.hyper_endpt=Vtimepoints[-1]-capacitive_artifact_points #make sure capacitive artifact excluded
            self.hyperstartpt=self.hyper_endpt-int(self.ss_dur/self.dt)
            self.hyperstart=self.hyperstartpt*self.dt
            self.Iaccess=inject[0]-np.mean(inject[timepoints[-2]:timepoints[-1]],axis=0) #use timepoints from current
            self.anal_params['hyperstart']=self.hyperstart
            self.params['Iaccess']=self.Iaccess[0]
            print('Access Monitor: V onset: {:g} and {:.5f} s, V offset: {:g} , inject {:.1f} pA'.format(Vtimepoints[-2],self.hyperstart, self.hyper_endpt,self.Iaccess[0]*1e12))
            return Vtimepoints
        
        if len(self.induction_time):
            induct_time=self.induction_time[0]
        else:
            induct_time=self.get_routine_start(self.pre_trace)[0] #count psptime from beginning of baseline?
        for headstage in self.Vm_psp.keys():  #intersection of headstages requested and those available
            trace_num=0 #trace_num indexes the measurement array
            for r in self.Vm_psp[headstage].keys():
                self.dt=self.get_dt(self.psp_dict[r])
                routine_start,routine_name=self.get_routine_start(r)
                #determine time and value of hyperpolarizing pulses - two pulses detected
                timepoints,inject=self.find_inject(r)
                Vtimepoints=access_monitor(timepoints,inject, self.dt, r)
                if isinstance(self.PSPstart,dict): #if PSPstart extracted from notebook file
                    self.PSPstartpt=int(self.PSPstart[routine_name]/self.dt)                    
                else: #if no notebook file
                    self.PSPstartpt=int(self.PSPstart/self.dt) 
                if self.PSPstartpt>np.shape(self.data['Data'][r].__array__())[0]:
                    print('********* PSP start time', self.PSPstart,' is after end of trace.  Make start time earlier ***********')
                self.basestartpt=self.PSPstartpt-int(self.basestart/self.dt)  #basestart is duration prior to event 
                self.base_endpt=self.basestartpt+int(self.base_dur/self.dt)
                inj2ms=inject[timepoints[1]]-inject[timepoints[0]]
                for num in range(np.shape(self.data['Data'][r].__array__())[-1]):  #num indexes arrays - resets to zero for each routine
                    trace=self.data['Data'][r].__array__()[:,num] 
                    self.psp[headstage][trace_num],self.RMP[headstage][trace_num],peakpt=self.psp_detect(r,num,self.dt,self.PSPstartpt)
                    self.pspamp[headstage][trace_num]=self.psp[headstage][trace_num]-self.RMP[headstage][trace_num]
                    self.peaktime[headstage][trace_num]=(peakpt)*self.dt 
                    if len(self.sweep_time[routine_name]):
                        self.psptime[headstage][trace_num]=self.sweep_time[routine_name][num]-induct_time
                    else:
                        self.psptime[headstage][trace_num]=routine_start+num*self.params['Stim_interval'][routine_name]-induct_time
                    #### Potential Problem: use RMP from beginning of trace as RMP 1/5 sec later.
                    self.Raccess[headstage][trace_num]=(self.RMP[headstage][trace_num]-np.mean(trace[self.hyperstartpt:self.hyper_endpt]))/self.Iaccess[num]
                    if len(Vtimepoints)>=4:
                        self.anal_params['dV_2ms_time']=((Vtimepoints[0]+1)*self.dt, (Vtimepoints[1]+2)*self.dt) #add 1 point to account for deriv, 1 point for current to actually stop
                        self.dV_2ms[headstage][trace_num]=(trace[Vtimepoints[1]+2]-trace[Vtimepoints[0]+1]) #/inj2ms
                    trace_num+=1

    def psp_detect(self, r, num, dt, stim_start):
        trace=self.data['Data'][r].__array__()[:,num]
        peakpt=np.argmax(trace[ stim_start:])+ stim_start
        maxvm=np.mean(trace[peakpt-self.window:peakpt+self.window])
        if maxvm>bad_psp: #an AP occurred.  This might not be best criteria 
            maxvm=np.nan 
        basestartpt=stim_start-int(self.basestart/dt)  #basestart is duration prior to event 
        base_endpt=basestartpt+int(self.base_dur/dt)
        RMP=np.mean(trace[basestartpt:base_endpt])
        return maxvm, RMP, peakpt

    def IO_curve(self): #better to separate out colors and analysis, but would need to save more values
        from matplotlib import pyplot as plt
        fig,axes=plt.subplots(len(self.IO_psp),1)
        axes=fig.axes
        colors=plt.get_cmap('plasma')
        partial_scale=0.9 #avoid the very light 
        offset=.002 #in mV, to visualize multiple traces
        self.IOamp={headstage:np.zeros(len(self.IOrange)) for headstage in self.headstages} #peak value of psp
        for hd,headstage in enumerate(self.IO_psp.keys()):  #intersection of headstages requested and those available
            for r in self.IO_psp[headstage].keys():
                dt=self.get_dt(self.IO_dict[r])
                routine_start,routine_name=self.get_routine_start(r)
                IOtraces=np.shape(self.data['Data'][r].__array__())[-1]
                if len(self.IOrange) != IOtraces:
                    print('samples in IO range,',len(self.IOrange), ', does not match number of IO traces, ',IOtraces, ', truncating IOrange to', self.IOrange[0:IOtraces])
                    self.IOrange=self.IOrange[0:IOtraces]
                    self.IOamp[headstage]=np.zeros(IOtraces)
                if isinstance(self.PSPstart,dict):
                    stim_start=int(self.PSPstart[routine_name]/dt)
                else:
                    stim_start=int(0.20/dt) #FIXME: make parameter in ArgParser, 0.203?
                #### plotting part ####
                basestartpt=stim_start-int(self.basestart/dt)  #basestart is duration prior to event 
                base_endpt=basestartpt+int(self.base_dur/dt)
                colinc=(len(colors.colors)-1)/(IOtraces-1)
                for num in range(IOtraces):  #num indexes arrays - resets to zero for each routine
                    maxvm,RMP,peakpt=self.psp_detect(r,num,dt,stim_start)
                    self.IOamp[headstage][num]=maxvm-RMP
                    ##### Plotting part ####
                    trace=self.data['Data'][r].__array__()[:,num] 
                    time=np.arange(len(trace))*dt
                    color_index=int(num*colinc*partial_scale)
                    mycolor=colors.colors[color_index]
                    axes[hd].plot(time,trace+num*offset, color=mycolor, label=self.IOrange[num])
                    axes[hd].plot(time[peakpt],trace[peakpt]+num*offset,'r*')
                    for pt in [basestartpt,base_endpt]:
                        axes[hd].plot(time[pt],trace[pt]+num*offset,'k.')
            axes[hd].set_xlabel('Time (sec)')
            axes[hd].set_ylabel(headstage+' Vm (volts)')
            axes[0].set_title ('IO curve')
            axes[0].legend()  
     
    def normalize(self,baseline_time=None):
        from scipy import optimize
        self.meanpre={}
        self.Aopt={};self.Bopt={}; self.Bstd={}
        self.num_pre=np.shape(self.data['Data'][self.pre_trace].__array__())[-1]
        if self.baseline_time:
            pre_length=int((self.baseline_time*60)/self.Stim_interval) #convert from minutes to seconds to samples
            pre_start=max(0,self.num_pre-pre_length) #cannot be negative
        else:
            pre_start=0
        for headstage in self.pspamp.keys():
            pretrace=self.pspamp[headstage][pre_start:self.num_pre]
            pretime=self.psptime[headstage][pre_start:self.num_pre]
            self.meanpre[headstage]=pretrace[~np.isnan(pretrace)].mean()
            self.normpsp[headstage]=self.pspamp[headstage]/self.meanpre[headstage]
            #normRMP or normRaccess? 
            if len(pretrace) > 1:   
                popt,pcov=optimize.curve_fit(pu.line,pretime[~np.isnan(pretrace)],pretrace[~np.isnan(pretrace)]) #NOTE: fitting to pretrace means fitting to unnormalized data
                self.Aopt[headstage],self.Bopt[headstage]=popt
                Astd,self.Bstd[headstage]=np.sqrt(np.diag(pcov))
                if  np.abs(self.Bopt[headstage])-self.slope_std_factor*self.Bstd[headstage]>0:
                    print( "********************** WARNING ************** \n Baseline Slope of ", round(self.Bopt[headstage],7),"+/- 2*",round(self.Bstd[headstage],7), "sig diff than zero")
                else:
                    print( "**********************Baseline Seems to Be Fine********************: slope=", round(self.Bopt[headstage],7))
            else:
                self.Bopt[headstage]=self.Bstd[headstage]=self.Aopt[headstage]=np.nan

    def write_data(self):
        outfname=self.outputdir+self.experiment
        data_dict={'slope':self.Bopt, 'Intercept':self.Aopt,'slope_std':self.Bstd,'meanpre':self.meanpre,'max_latency':self.max_latency,'rheobase':self.rheobase,'psp_dt':self.dt, 'num_pre':self.num_pre} #single values
        tracedict={'amp':self.pspamp,'RMP':self.RMP, 'Raccess': self.Raccess,'normPSP': self.normpsp,'psptime':self.psptime,'peaktime':self.peaktime,'dV_2ms':self.dV_2ms} #arrays
        IV_IFsummary={'IV':self.IV_IF,'spikes':self.IV_IF_spikes,'dt':self.IV_dt} #analysis of IV_IF routines
        IOsummary={'amp':self.IOamp}
        if len(self.induction_list):
            induct_dict={'spikes':self.num_spikes,'induction':self.params['induction'], 'ISImin':self.min_interval, 'ISIother':self.other_interval, 'time_induct':self.params['time_to_induct'], 'spikes':self.spikes} #in GrpAvg, calculate mean spikes per routine?
            np.savez(outfname, trace=tracedict,params=self.params,data=data_dict,IV_IF=IV_IFsummary,anal_params=self.anal_params,IO=IOsummary,induct=induct_dict)
        else:
            np.savez(outfname, trace=tracedict,params=self.params,data=data_dict,IV_IF=IV_IFsummary,anal_params=self.anal_params,IO=IOsummary)

if __name__=='__main__':
    #ARGS='250704_1 -headstages H2 -celltype D1-SPN -decay .003 -base_time 5'
    try:
        commandline = ARGS.split() 
        do_exit = False
    except NameError: 
        commandline = sys.argv[1:]
        do_exit = True

    params=argp.ArgParserPatch(commandline,do_exit)
    print(commandline, params)
    exp=PatchAnal(params)
    exp.read_datafile()
    exp.read_notebook()
    exp.time() #May not be needed if Notebook file saved, and if time is in notebook for every routine
    exp.meta_data()
    exp.init_arrays()
    exp.IV_IF_anal()
    exp.calc_psp()
    exp.normalize() #Specify baseline_time in minutes if not entire baseline data
    if len(exp.induction_list):
        stim_time=exp.analyze_theta()
    exp.IO_curve()
    exp.write_data()
    if exp.graphs:
        pu.IVIF_plot(exp)
        pu.IVIF_measures(exp)
        pu.trace_plot(exp) #set ss_color='k' if you want black, etc.  set ss_color=None to turn off lines
        pu.summary_plot(exp)
        pu.IO_plot(exp)
        if len(exp.induction_list):
            pu.induction_plot(exp,stim_time)
    print('time to induction=', round(exp.params['time_to_induct']/60,2),'min')



