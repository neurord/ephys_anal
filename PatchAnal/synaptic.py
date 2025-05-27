import numpy as np
import sys
import argparse
import pandas as pd
from spike_utilities import find_notebook_file
from PatchAnalPlots import rows_columns, plot_hist, plot_peaks, plot_pairs

HEADSTAGE_I={'S1': 'H1', 'S2': 'H2'}

def ArgParserSyn(commandline,do_exit):
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', type=str, help = 'give name of experiment')
    parser.add_argument('sweeps', type=int, nargs="+", help='start and end sweep numbers') 
    parser.add_argument('-datadir',type = str, default = "\\\\iowa.uiowa.edu\\shared\\engineering\\home\\rrajn\\windowsdata\\Desktop\\Recordings\\",   help="Directory containing raw data")
    parser.add_argument('-headstage', type=str, nargs="+", help='which headstage to analyze (default=both)', default=['H1','H2']) #FIXME - make this required for PatchAna?
    parser.add_argument("-celltype", type=str, nargs="+", choices=["D1-SPN", "FSI", "D2-SPN",'ChI','other'],help="D1-SPN, D2-SPN, FSI or other cell type, one for each headstage") #FIXME - make this required for PatchAnal?
    parser.add_argument("-bins", type=int,help="Number of bins for probability distributions", default=40) 
    parser.add_argument("-use_h5", type=str,help="use h5 file to extract events", default='y') 
    parser.add_argument("-win", type=int,help="window for filtering and finding peak value for halfwidth", default=5)
    parser.add_argument("-base", type=float,help="width for calculating baseline", default=.004)
    parser.add_argument("-decay", type=float,help="duration to search for decay time", default=.03)
    parser.add_argument("-filter", type=int,help="cutoff frequency of 4 pole butterworth lowpass, set to 0 to avoid filtering", default=1000)
    parser.add_argument("-height", type=float,help="minimum amplitude of PSC peaks", default=13e-12)
    parser.add_argument("-wlen", type=float,help="distance to search for PSC peaks", default=0.06)
    parser.add_argument("-IPSC_routine",type=str,default='R2')
 
    try:
        args = parser.parse_args(commandline) # maps arguments (commandline) to choices, and checks for validity of choices.
    #if arguments are mapped incorrectly, python wants to exit, but the next line says "don't", instead check whether we are in python (do_exit=False) then don't exit, just give us a warning
    except SystemExit:
        if do_exit:
            raise # raise the exception above (SystemExit) b/c none specified here
        else:
            raise ValueError('invalid ARGS')
    return args

class SynAnal():
    def __init__(self,params): #change to PopSpike(args)
        self.datadir = params.datadir
        #self.outputdir = params.outputdir
        self.exper=params.exper
        if len(params.headstage)==len(params.celltype): 
            self.celltypes={h:params.celltype[i] for i,h in enumerate(params.headstage)}
        else:
            print('*********** Need to specify celltype for each valid headstage ***********')
            exit()
        self.params={'exper':self.exper}
        self.bins=params.bins
        self.win=params.win
        self.baseline_width=params.base
        self.filt_cutoff=params.filter
        self.height=params.height
        self.wlen=params.wlen
        self.decay_time=params.decay
        self.dt=0
        if params.use_h5=='y' or params.use_h5=='Y':
            self.h5file=params.datadir+self.exper+'.h5'
        self.IPSC_routine_num=params.IPSC_routine
        
    def read_notebook(self): #extract routine and sweep time from notebook file, if it was saved
        def parse_lines(all_lines):
            event_start={};event_end={};column_names=[]
            for i,line in enumerate(all_lines):
                if line.startswith('Mini Results'): #start of synaptic data for one headstage                    
                    lin=line.strip('\n')
                    event_start={lin:i}                    
                    while event_start[lin]<len(all_lines): #keep reading until the first numeric line
                        entry=all_lines[event_start[lin]].split(',')[0].strip(' ')
                        if entry.isnumeric():
                            event_end={lin:event_start[lin]}
                            break
                        else:
                            event_start[lin]+=1
                    for j,ln in enumerate(all_lines[event_start[lin]:]): #Now, keep reading till the end of the data for that headstage
                        if ln.startswith('\n'):
                            event_end[lin]+=j
                            break
                    channel=line.split('_')[-1].strip('\n')
                    headstage=HEADSTAGE_I[channel]
                    celltype=self.celltypes[headstage]
                elif line.startswith('"Entry"'):
                    full_line=''.join([all_lines[j].strip('\n') for j in range(i,event_start[lin])])
                    columns=full_line.split(',')
                    column_names=[c.strip(' ').strip('"') for c in columns]
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
            return event_start,event_end, column_names

        unit_dict={'m':1e-3,'u':1e-6,'n':1e-9,'p':1e-12,'f':1e-15, 'a':1e-18}
        def value_units(line, columns):
            temp_dict={}
            for i,entry in enumerate(line):
                parts=entry.strip().split()
                value=parts[0]
                #convert to either float or int
                if '.' in value:
                    value=float(value)
                elif value.isnumeric():
                    value=int(value) 
                elif value=='NaN':
                    value=np.nan
                #extract units, if there are units associated with the entry
                if len(parts)>1:
                    units=parts[-1]
                    #convert to SI units
                    if units == 's' or units == 'xSD' or units=='px1e-12':
                        value=value #essentially do nothing
                    else:
                        value=value*unit_dict[units[0]]
                temp_dict[columns[i]]=value
            return temp_dict

        import glob
        import os
        nfname=self.datadir+'Notebook_20'+'*.txt'  #find all Notebook files with correct date
        files=glob.glob(nfname)
        if len(files):
            findex,file_found=find_notebook_file(files,self.exper)
            file_found=True
        else:
            file_found=False
        if file_found: 
            print('Using notebook file:',files[findex])
            self.params['file']=os.path.basename(files[findex]).split('.')[0]  
            with open(files[findex],'r') as myfile: 
                all_lines=myfile.readlines()
            event_start,event_end,column_names=parse_lines(all_lines)
            self.dfset={}
            for key in event_start.keys(): #will need event_end when there are two headstages
                PSC_list=[]
                for line in range(event_start[key],event_end[key]):
                    data=all_lines[line].split(',')                  
                    PSC_list.append(value_units(data, column_names))
                self.dfset[key] = pd.DataFrame.from_dict(PSC_list)
                self.dfset[key]['ISI_orig']=self.dfset[key]['Interevent Interval (s)'] 
                self.dfset[key]['Interevent Interval (s)']=np.insert(np.diff(self.dfset[key]["Absolute Event Time (s)"]),0,np.nan)
        else:
            print('0 or multiple Notebook files found using pattern:',nfname)
        return

    def analyze(self,measures,dfset):
        self.ecdf={k:{} for k in dfset.keys()}
        #self.bin_edges={k:{} for k in self.dfset.keys()}
        results={k:{} for k in dfset.keys()}
        self.hist_measures={'IEI':'Interevent Interval (s)',
                  'amp':'Event Amplitude (A)'}
        self.units={'(A)':'1e12','(A*s)':'1e12'} #do not list units conversions already in seconds
        for key in dfset.keys():
            for m,nm in measures.items():
                if nm.split()[-1] in self.units.keys():
                    dfset[key][nm]=dfset[key][nm]*float(self.units[nm.split()[-1]]) #convert from A to pA
                results[key][m+'_mean']=dfset[key].mean()[nm]
                results[key][m+'_stdev']=dfset[key].std()[nm]
            results[key]['Freq_mean']=np.nanmean(1/dfset[key]['Interevent Interval (s)']) #FIxME: 
            results[key]['Freq_std']=np.nanstd(1/dfset[key]['Interevent Interval (s)']) #FIxME: 
            #calculate histogram for export - may not be needed
            for meas,long_name in self.hist_measures.items():
                if meas=='IEI':
                   self.ecdf[key][meas]=np.diff(dfset[key]["Absolute Event Time (s)"]) #eliminate nans at beginning of each sweep
                else:
                    self.ecdf[key][meas]=dfset[key][long_name].values
                #minv=np.percentile(df[long_name].dropna().values,1)
                #maxv=np.percentile(df[long_name].dropna().values,99)
                #self.hist[meas],self.bin_edges[meas]=np.histogram(df[long_name].values,range=(minv,maxv),bins=self.bins)
            return results

    def read_h5(self,key=None,headstage=None):
        import h5py as h5
        if key:
            self.psc_name=key.split('Results_')[-1] 
        else:
            self.psc_name=self.IPSC_routine_num+'_'+list(HEADSTAGE_I.keys())[list(HEADSTAGE_I.values()).index(headstage)] +'_IPSC'           
        self.data = h5.File(self.h5file,"r")
        self.routines=list(self.data['Data'].keys())
        self.trace_name=[r for r in self.routines if r.startswith(self.psc_name)][0]
        self.get_dt(self.trace_name)
 
    def halfwidth(self): #assumes outward (upward) currents.  FIXME: add param - if inward, invert trace_subset
        import scipy.signal        

        for ky in self.dfset.keys():
            if not self.dt:
                self.read_h5(ky)
            halfw=[];decay=[]
            trace_subset=self.data['Data'][self.trace_name].__array__()[:,self.dfset[ky]["Sweep Number"].min()-1:self.dfset[ky]["Sweep Number"].max()]
            fig,axes=self.plot_traces(trace_subset,ky=ky,title='dPatch Events')
            b,a=scipy.signal.butter(4,1000,btype='lowpass',fs=(1./self.dt)) #1 khz 4 pole filter
            trace_filtered=scipy.signal.filtfilt(b,a, trace_subset, axis=0)
            axes[0].plot(np.arange(len(trace_subset))*self.dt,trace_filtered[:,0]) 
            mean_rise_dur=3*self.dfset[ky]['10-90% Rise Time (s)'].mean() #FIXME: factor of 2 or 3 is arbitrary
            for ind,tm in enumerate(self.dfset[ky]["Event Time (s)"]):
                point=int(tm/self.dt) #peak point from dpatch
                trace_num=self.dfset[ky]["Sweep Number"][ind]-self.dfset[ky]["Sweep Number"].min()
                #decay_dur is duration to search for half max during decay phase, cannot go beyond next event or end of recording
                if ind<len(self.dfset[ky]["Event Time (s)"])-1:
                    decay_end=int((self.dfset[ky]["Absolute Event Time (s)"][ind+1]-self.dfset[ky]["Absolute Event Time (s)"][ind])/self.dt) 
                else:
                    decay_end=np.shape(trace_subset)[0]-peakpt
                decay_dur=min(int(self.decay_time/self.dt),decay_end) #number of points to search for half decay after peakpt
                if self.dfset[ky]['10-90% Rise Time (s)'][ind]<=0:
                    self.dfset[ky]['10-90% Rise Time (s)'][ind]=0.2e-3
                rise_dur=int(mean_rise_dur/self.dt)#int(4*self.dfset[ky]['10-90% Rise Time (s)'][ind]/self.dt)  #duration to search for half max during rising phase 
                peakpt=trace_filtered[point-rise_dur:point+decay_dur,trace_num].argmax()+point-rise_dur  #better estimate of peak location?
                #base_end=np.diff(trace_filtered[point-rise_dur:point,trace_num]).argmax()+point-rise_dur
                base_end=point-rise_dur #point is usually earlier than peakpt
                base_start=max(base_end-int(self.baseline_width/self.dt),0)
                baseline=trace_subset[base_start:base_end,trace_num].mean() 
                peakval=trace_filtered[peakpt,trace_num] #get peak from filtered trace
                halfmax=(peakval-baseline)/2+baseline #half amplitude of PSC
                halfup_points=np.where(trace_filtered[peakpt-rise_dur:peakpt,trace_num]<halfmax)
                if len(halfup_points[0]):
                    halfup=np.max(halfup_points)+peakpt-rise_dur#halfway point during rising phase
                else:
                    halfup=peakpt-rise_dur  
                halfdown_points=np.where(trace_filtered[peakpt:peakpt+decay_dur,trace_num]<halfmax) #halfway point during decay phase
                if len(halfdown_points[0]):
                    halfdown=np.min(halfdown_points)+peakpt  #halfway point during rising phase
                else:
                    halfdown=peakpt+decay_dur
                halfw.append((halfdown-halfup)*self.dt)
                decay_val=baseline+0.37*(peakval-baseline)
                decay_points=np.where(trace_filtered[peakpt:peakpt+decay_dur,trace_num]<decay_val)
                if len(decay_points[0]):
                    decay_pt=np.min(decay_points)+peakpt
                else:
                    decay_pt=decay_dur+peakpt
                decay.append((decay_pt-peakpt)*self.dt)
                #plot peak
                axes[trace_num].plot(peakpt*self.dt,peakval,'r*')
                #plot baseline
                xval=np.arange(base_start,base_end)*self.dt
                axes[trace_num].plot(xval,[baseline]*len(xval),'cyan')
                #plot half width
                xval=np.arange(halfup,halfdown)*self.dt
                #axes[trace_num].plot(xval,trace_filtered[halfup:halfdown,trace_num],'cyan')
                axes[trace_num].plot(xval,[halfmax]*len(xval),'m')
                #visualize the trace - FIXME: 50 ms before and 100 ms after peak is arbitrary
                axes[trace_num].set_xlim([peakpt*self.dt-.050,peakpt*self.dt+0.100])
                axes[trace_num].plot(decay_pt*self.dt,trace_filtered[decay_pt,trace_num],color='purple',marker='*',ms=8)
                axes[trace_num].plot(peakpt+decay_dur,trace_filtered[int(self.decay_time/self.dt)+peakpt,trace_num],color='black',marker='|',ms=15)
            axes[0].set_xlim([0,len(trace_subset)*self.dt])
            self.dfset[ky]['width (s)']=halfw
            self.dfset[ky]['Event Decay Tau (s)']=halfw
        #
    def detect_events(self, trace_nums,width=5e-3,polarity='plus'): #width used to define distance between peaks and minimum width of peaks
        from scipy.signal import find_peaks, butter, filtfilt
        import h5py as h5
        self.dfset2={}
        for ky in self.celltypes.keys():
            if not self.dt:
                self.read_h5(headstage=ky)
            traces=self.data['Data'][self.trace_name].__array__()[:,trace_nums[0]-1:trace_nums[1]]
            peaks=[]
            if self.filt_cutoff>0:
                b,a=butter(4,self.filt_cutoff,btype='lowpass',fs=(1./self.dt)) # 4 pole filter
                trace_filt=filtfilt(b,a, traces, axis=0)
            else:
                trace_filt=traces
            for trace_num in range(np.shape(traces)[1]):
                if polarity=='minus':
                    mult=-1
                else:
                    mult=1
                #trace_filt.append(uniform_filter1d(mult*traces[:,trace_num],self.win))
                #try wlen=100 ms?
                peaks.append(find_peaks(mult*trace_filt[:,trace_num],prominence=self.height,width=int(width/self.dt),distance=int(width/self.dt),wlen=int(self.wlen/self.dt)))
            fig,axes=self.plot_traces(traces,peaks=peaks,trace_filt=trace_filt)
            event_amp=[[traces[jj][trace_num] for jj in peaks[trace_num][0]] for trace_num in range(len(peaks))]
            event_time=[peaks[trace_num][0]*self.dt for trace_num in range(len(peaks))]
            plot_peaks(fig,axes,event_time, event_amp,'k',symbol='o')
            cont=input(str(np.sum([len(ps[0]) for ps in peaks]))+' peaks found, continue (y/n, if n, program will exit)?') 
            if cont=='n' or cont == 'N':
                sys.exit()
            PSC_list=self.check_peaks(trace_filt,peaks) #use trace_filt for better estimation of amplitude and baseline
            plot_peaks(fig,axes,self.peak_data['event_time'].values(),self.peak_data['amp'].values(),'r',symbol='*') 
            fig.suptitle('black dot=detected peaks, red * = accepted peaks')
            self.dfset2[ky]=pd.DataFrame.from_dict(PSC_list) 
            self.dfset2[ky]['Interevent Interval (s)']=np.insert(np.diff(self.dfset2[ky]['Absolute Event Time (s)']),0,np.nan)

    def check_peaks(self,traces,peaks,surround=0.05): #50 ms before and after for visualization is arbitrary
        from matplotlib import pyplot as plt
        #1. fix baseline:
        self.new_measures={'event_time', 'amp'}
        self.peak_data={m:{x:[] for x in range(len(peaks))} for m in self.new_measures}
        PSC_list=[]
        for sn,peak_set in enumerate(peaks):
            print('###########', len(peak_set[0]), 'peaks detected in sweep ', sn)
            fig,ax=plt.subplots(1,1,figsize=(4,2))
            time=np.arange(len(traces))*self.dt
            ax.plot(time,traces[:,sn])
            ax.set_title('trace '+str(sn))
            for ii,pt in enumerate(peak_set[0]): #peak_set[0] is list of points.  peak_set[1] is dict of peak properties
                #adjust xaxis to show single event
                min_tm=max(pt*self.dt-surround,0)
                max_tm=min(pt*self.dt+surround,time[-1])
                ax.set_xlim([min_tm,max_tm]) #just show a single peak.  Color code so it stands out
                peak_val=traces[pt,sn]
                #how to define start of PSC - steepest slope from derivative?, min point with 5 ms of left halfwidth
                #how to define decay time?  reach within 10% of baseline? time from peak to half decay - well defined
                left_half_height=int(peak_set[1]['left_ips'][ii])
                start=max(0,left_half_height-int(self.baseline_width/self.dt)) #go backward base_width from left width point to search for PSC start
                #define psc start as steepest derivative within 2 ms prior to left half height
                base_end=peak_set[1]['left_bases'][ii]
                #define psc start as minimum point within 2 ms prior to left half height
                #print('base end deriv=',base_end,round(base_end*self.dt,4),end="")
                base_end=np.argmin(traces[start:left_half_height,sn])+start #find minimum
                #print(', *base end min=',base_end,round(base_end*self.dt,4),'base end "bases"=',round(peak_set[1]['left_bases'][ii]*self.dt,4),end="")
                #define baseline starting point as base_width ms prior to base_end/psc_start
                base_start=max(base_end-int(self.baseline_width/self.dt),0)
                baseline=np.mean(traces[base_start:base_end,sn]) #measure baseline prior to location of minimum used for amplitude
                #print(' baseline=',baseline,'baseline from find_peaks=',max(traces[peak_set[1]['left_bases'][ii],sn],traces[peak_set[1]['right_bases'][ii],sn]))
                amp=peak_val-baseline #alternative amplitude: peak to baseline
                if len(np.where(traces[base_end:pt,sn]<amp*0.1+baseline)[0]):
                    rise10=np.max(np.where(traces[base_end:pt,sn]<amp*0.1+baseline))+base_end
                else:
                    rise10=base_end #search from the minimum to peak for 10% rise
                if len(np.where(traces[base_end:pt+1,sn]>amp*0.9+baseline)[0]):
                    rise90=np.min(np.where(traces[base_end:pt+1,sn]>amp*0.9+baseline))+base_end
                else:
                    rise90=pt #search from the minimum to peak for 90% rise
                decay_val=baseline+0.37*amp
                decay_points=np.where(traces[pt:int(self.decay_time/self.dt)+pt,sn]<decay_val) #FIXME: limit decay search to time of next event
                if len(decay_points[0]):
                    decay_pt=np.min(decay_points)+pt
                else:
                    decay_pt=int(self.decay_time/self.dt)+pt
                decay=(decay_pt-pt)*self.dt
                ######### plot width, peak,baseline on graph ###################
                #plot width at half height
                xvals=time[int(peak_set[1]['left_ips'][ii]):int(peak_set[1]['right_ips'][ii])]
                ax.plot(xvals,[peak_set[1]['width_heights'][ii]]*len(xvals),'m')
                #plot peak
                ax.plot(time[pt],peak_val,'ko')
                ax.plot(time[base_start:base_end],[baseline]*len(np.arange(base_start,base_end)),'red')
                ax.plot(time[rise10:rise90],traces[rise10:rise90,sn],'cyan')
                ax.plot(time[decay_pt],traces[decay_pt,sn],color='purple',marker='*',ms=8)
                ax.plot(time[int(self.decay_time/self.dt)+pt],traces[int(self.decay_time/self.dt)+pt,sn],color='black',marker='|',ms=15)
                keep=input('keep this PSC from sweep '+str(sn)+' at '+str(round(pt*self.dt,4))+'s (y/n)? ')
                #keep='y'
                #for art in list(ax.lines[1:]):art.remove() #remove the color coding and peak point?
                if keep=='y' or keep == 'Y':
                    mydict={'Event Amplitude (A)':amp,'Absolute Event Time (s)':pt*self.dt+sn*time[-1],
                            '10-90% Rise Time (s)':(rise90-rise10)*self.dt,'width (s)':peaks[sn][1]['widths'][ii]*self.dt,
                            'prom (A)':peaks[sn][1]['prominences'][ii],'Event Decay Tau (s)':decay}
                    PSC_list.append(mydict)
                    self.peak_data['amp'][sn].append(peak_val)
                    self.peak_data['event_time'][sn].append(pt*self.dt)
                    #self.peak_data['rise'][i].append((rise90-rise10)*self.dt)
                    #self.peak_data['width'][i].append(peaks[i][1]['widths'][ii]*self.dt)
                    #self.peak_data['prom'][i].append(peaks[0][1]['prominences'][ii])
        return PSC_list

    #keep in class because this plotting function has several class variables
    def plot_traces(self,trace_subset,ky=None,peaks=None,trace_filt=None,title=''):
        from matplotlib import pyplot as plt
        time=np.arange(len(trace_subset))*self.dt
        fig,axes=plt.subplots(nrows=np.shape(trace_subset)[1],ncols=1,sharex=True,figsize=(12,6))#(22,13))
        for trace_num in range(np.shape(trace_subset)[1]):
            axes[trace_num].plot(time,trace_subset[:,trace_num])
            axes[trace_num].set_ylabel('trace '+str(trace_num))
            if trace_filt is not None:
                axes[trace_num].plot(time,trace_filt[:,trace_num],color='y')
            if peaks:
                for i,(pt,width) in enumerate(zip(peaks[trace_num][0],peaks[trace_num][1]['widths'])): 
                    #plot width of event
                    left=peaks[trace_num][1]['left_ips'][i]
                    right=peaks[trace_num][1]['right_ips'][i]
                    xval=np.arange(left,right)*self.dt
                    yval=[peaks[trace_num][1]['width_heights'][i]]*len(xval)
                    axes[trace_num].plot(xval,yval,'magenta')
            if ky:
                trace_events=self.dfset[ky].loc[self.dfset[ky]['Sweep Number']==trace_num+self.dfset[ky]["Sweep Number"].min()]
                startevent=[int((trace_events['Event Time (s)'][p]-1.5*trace_events['10-90% Rise Time (s)'][p])/self.dt) for p in trace_events.index]
                base_dur=int(0.002/self.dt) #FIXME: make this number a parameter.  1.5 in previous line is also arbitrary - prob not sufficient
                baseline=[np.mean(trace_subset[rp-base_dur:rp,trace_num]) for rp in startevent]
                axes[trace_num].scatter(trace_events['Event Time (s)'],trace_events['Event Amplitude (A)']+baseline,marker='o',color='k')
        return fig,axes

    def get_dt(self,routine_name): #use this to find dt for any routine 
        routine_list=[]
        for row in self.data['DataWaveNotes'].__array__(): #units, number of points, min and max values, dt for traces
            routine_list.append(row[0].decode('UTF8'))
        routine_num=routine_list.index(routine_name)
        #R1,R2, etc = number of routine
        #S1: channel 1 current,S2: channel 1 voltage;S3-channel2 current;S4-channel 2 current
        text=self.data['DataWaveNotes'].__array__()[routine_num][2].decode('UTF8')
        self.dt=float(text.split('XINC=')[-1].split(';')[0]) #needs to be identifical for both headstages and all psp routines

    def save_data(self,results,dfset):
        self.params['SxDate']=self.SurgeryDate
        self.params['ID']=self.ID
        self.params['region']=self.region
        self.params['drug']=self.drug
        if 'exper' in self.params:
            outfname=self.datadir+self.params['exper']
        else:
            outfname=self.datadir+self.params['file']
        for key in dfset.keys():
            if key in HEADSTAGE_I.values():
                h=key
            else:
                h=HEADSTAGE_I[key.split('_')[-1]]
            self.params['celltype']=self.celltypes[h]
            #FIXME: find_peaks file will over-write the npz file using dPatch peaks
            np.savez(outfname+h, params=self.params,data=results[key],cdf=self.ecdf[key])


if __name__=='__main__':
    #uncomment (and edit) to run in debugger
    #ARGS='250501_4 6 8 -headstage H2 -celltype D1-SPN -bins 50 -height 15e-12 -wlen 0.08 -datadir IPSCs/'
    try:
        commandline = ARGS.split() 
        do_exit = False
    except NameError: 
        commandline = sys.argv[1:]
        do_exit = True

    params=ArgParserSyn(commandline,do_exit)
    print(commandline, params)
    exp=SynAnal(params)
    exp.read_notebook()
    measures={'IEI':'Interevent Interval (s)', 
            'rise':'10-90% Rise Time (s)', 
            'decay':'Event Decay Tau (s)',
            'amp':'Event Amplitude (A)', 'width':'width (s)'}
    #use events detected by dPatch, with values written to notebook; need to measure halfwidth
    '''    exp.halfwidth()
        for key,df in exp.dfset.items():
            print(key,'           CORRELATIONS for peaks found with dPatch\n',df.corr(),'\n')

        results=exp.analyze(measures,exp.dfset,)
        exp.save_data(results,exp.dfset)
    '''
    #peak detection using scipy's find_peaks
    exp.detect_events(params.sweeps,width=5e-3) #image size, height of events, width used for detection
    results2=exp.analyze(measures, exp.dfset2)
    exp.save_data(results2,exp.dfset2)
    
    for key,df in exp.dfset2.items():
        print(key,'           CORRELATIONS for peaks found with scipy.find_peaks\n',df.corr(),'\n')
    fighist=plot_hist(measures.values(),exp.dfset2,exp)
    pairs=[('Event Amplitude (A)','Event Decay Tau (s)'),('Event Decay Tau (s)','10-90% Rise Time (s)')]
    fig=plot_pairs(exp.dfset2,pairs)


    
