import numpy as np
import argparse
import pickle

def parse_args(commandline,do_exit,flag):
    #Descriptive variables used in subsequent program to select similar data to average:
    parser = argparse.ArgumentParser()
    exp=['exper','--exper']
    sex=['sex','--sex']
    age=['age','--age']
    drug=['drug','--drug']
    theta=['theta','--theta']
    region=['region','--region']
    estratiol=['estradiol','--estradiol']
    parser.add_argument(exp[flag], type=str, help = 'give exp name, for example: 031011_5Or')
    parser.add_argument('--no-graphs', '-g', dest='graphs', default=True, action='store_false') # -g optional and not position-defined (its absence OK too)
    parser.add_argument(sex[flag], type=str, choices=["M","F", "Fe"],help="male M or female F or Fe")
    parser.add_argument(age[flag], type=int, help="animal age in days")
    parser.add_argument(drug[flag], type=str, help="what drugs were in the ACFS")
    parser.add_argument(theta[flag], type=float, help="theta frequency (e.g. 5, 8, 10.5), enter 0 for no stim ctrl")   
    parser.add_argument(region[flag], type=str, choices=["DM", "DL"],help="dorsomedial: DM or dorsolateral: DL")
    parser.add_argument(estradiol[flag], type=float, help="estradiol value, will default to -1 if value not specified", default=-1)
    
    # parser.add_argument('-sample_times', '--list', default = [30,60,90,100], action='append', help = "Enter to modify sample_times (default = [30,60,90,100]),", required = True) 
    parser.add_argument("--samp_time", nargs="+",default=[30,60,90,100],)   
    parser.add_argument("--sepvarlist", nargs="+",default=['sex','theta'],help='list of separation variables for grouping data')
    parser.add_argument('-decay', type = float,  default =1.7, help="artifact decay time (default: 1.7 msec)")
    parser.add_argument('-artthresh',  type = float, default = 1.5, help="artifact threshold (default: 1.5 mV")
    parser.add_argument('-FVwidth',  type = float, default = 1, help="Fiber volley width (default: 1 msec")
    parser.add_argument('-samp_win',  type = int, default =5, help="Average over this many minutes for time samples (default: 5)")
    parser.add_argument('-file_end',  type = str, default = ".Lvm", help="filename ending (default: .Lvm)")
    parser.add_argument('-first_peak_end_frac',  type = float, default = 0.625, help="first_peak_end_fraction (default: 0.625)")
    parser.add_argument('-art_win', type = float, default = 0.5, help="artifact window (default: 0.5)")
    parser.add_argument('-base_min', type = int, default = 15, help="number of minutes for stable baseline (default: 15)")
    parser.add_argument('-noisethres', type = float, default = 1.5, help="Noise threshold (default: 1.5)")
    parser.add_argument('-slope_std', type = int, default = 2, help="Enter to modify slope_std_factor (default: 2)")
    parser.add_argument('-datadir',type = str, default = ".\\",   help="Directory containing raw data")
    parser.add_argument('-outputdir',  type = str, default = ".\\", help="Directory for pickle files")
    if flag:
        parser.add_argument('--maxage', type=int, help="maximum age of animals")
    keydict={"sex":["M","F","Fe"],"region":["DM", "DL"]}

    with open('choicedict.txt', "wb") as f: #open the choice dictionary, allowing it to be rewritten
        pickle.dump(keydict, f) # same the choice dictionary, for use in analyzing groups of experiments

        try:
            args = parser.parse_args(commandline) # maps arguments (commandline) to choices, and checks for validity of choices.
        #if arguments are mapped incorrectly, python wants to exit, but the next line says "don't", instead check whether we are in python (do_exit=False) then don't exit, just give us a warning
        except SystemExit:
            if do_exit:
                raise # raise the exception above (SystemExit) b/c none specified here
            else:
                raise ValueError('invalid ARGS')
        
    #print experiment characteristics to double check unconstrained entries
    print ("exper={}".format(args.exper))
    print ("sex={}".format(args.sex))
    print ("age={}".format(args.age))
    print ("drug={}".format(args.drug))
    print ("theta={}".format(args.theta))
    print ("region={}".format(args.region))
    return args

def line(x,A,B):
    return B*x+A

def plot_peaks(exper,time,Vm_traces,peaktime,pospeaktime,peak,pospeak,popspikestart,base,goodtraces,baseline_start,FVwidth):
    traces_per_panel=40     #controls how many traces shown on each panel
    peak_decay=0.016        #controls how much of the trace is plotted after the peak
    spread=0.1              #controls how far apart to spread the traces 
    dt=time[1]-time[0]
    panels=int(len(peaktime)/traces_per_panel)
    if len(peaktime)%traces_per_panel:
            panels=panels+1
    fig,axes=pyplot.subplots(1,panels)
    fig.canvas.set_window_title('Good traces for '+exper)
    for panel in range(0,panels):
        axes[panel].clear()
        trace1=panel*traces_per_panel
        tracen=min(len(peaktime),panel*traces_per_panel+traces_per_panel)
        #fig.canvas.set_window_title('Traces '+str(trace1)+' to '+str(tracen)+' of '+exper)
        for index in range(trace1,tracen):
            if ~np.isnan(base[index]):
                st=baseline_start[index]
                #st=0 ##### uncomment to display entire trace
                end=int((peaktime[index]+peak_decay)/dt)
                #end=len(Vm_traces[goodtraces[index]]) #####uncomment to display entire trace
                offset=0.001-baseline_start[index]*dt  #display 1msec prior to baseline_start
                axes[panel].plot(time[st:end]+offset,Vm_traces[goodtraces[index],st:end]+index*spread,label=goodtraces[index])
                axes[panel].plot(peaktime[index]+offset,peak[index]+index*spread,'k*')
                if np.isnan(pospeak[index]):
                        #print "nan detected", index, pospeaktime[index]+offset,int(pospeaktime[index]/dt)
                        axes[panel].plot(pospeaktime[index]+offset,Vm_traces[goodtraces[index],int(pospeaktime[index]/dt)]+index*spread,'mD')
                elif base[index]>pospeak[index]:
                        axes[panel].plot(pospeaktime[index]+offset,pospeak[index]+index*spread,'r*')
                else:
                        axes[panel].plot(pospeaktime[index]+offset,pospeak[index]+index*spread,'ro')
                axes[panel].plot(dt*baseline_start[index]+offset,base[index]+index*spread,'bo')
                axes[panel].plot(popspikestart[index]+FVwidth+offset,Vm_traces[goodtraces[index],int(dt*popspikestart[index])]+index*spread,'b|')
                axes[panel].plot(popspikestart[index]+offset,Vm_traces[goodtraces[index],int(dt*popspikestart[index])]+index*spread,'k|')
        axes[panel].set_xlabel('Time (sec)')
        axes[panel].legend(fontsize=7, loc='right')
    fig.suptitle('black star=popspike,red o,* = pospeak, blue o = baseline, |=ps start/FV start')
    axes[0].set_ylabel('Vm (mV)')
    fig.canvas.draw()
            #text=raw_input('next? (y/n)')
    return fig,axes

def plot_summary(popspikeminutes,popspikenorm,baselineminutes,popspike_timesamples,start_samples,Aopt,Bopt,fig,axes):
    fig,axes=pyplot.subplots()
    fig.canvas.set_window_title('Summary ')
    axes.plot(popspikeminutes,popspikenorm,'b.')
    axes.plot(popspikeminutes[0:baselineminutes],line(popspikeminutes[0:baselineminutes],Aopt,Bopt),'r')
    start=np.insert(start_samples,0,0)
    time_summary=popspikeminutes[start]
    print ("summary at", len(start), "time points \n",np.column_stack((time_summary,popspike_timesamples[0:len(start)])))
    axes.plot(time_summary,popspike_timesamples[0:len(start)],'ko')
    return
        
def read_labview(filename):
    with open(filename,'r') as f:
            for line in f:
                    if "Samples" in line:
                            break
    timepoints_per_trace=int(line.split()[1])
    data = np.loadtxt(filename, skiprows=24) 
    #first column of datafile is time, second column is Vm
    wholetime=data[:,0]
    tempVm=data[:,1]
    dt=wholetime[1]-wholetime[0]
    xtime=wholetime[0:timepoints_per_trace]
    #reshape the data to put each trace in separate column
    datalength=np.shape(tempVm)[0]
    numtraces=int(datalength/timepoints_per_trace)
    Vm_traces=np.reshape(tempVm, (numtraces, -1))
    #extract the start time for each trace
    tracetime=np.zeros(numtraces)
    for i,j in enumerate(range(0, datalength, timepoints_per_trace)):
            tracetime[i]=wholetime[j]
    return Vm_traces,tracetime,dt,xtime
