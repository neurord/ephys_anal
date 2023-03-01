import numpy as np 
import glob
import pickle
import sys
import PopSpikeAnalClass as psc
import os
import time
import datetime
from argparse import Namespace
import pop_spike_utilities as psu
from PopSpikeReExtractPeaks import get_params

ms_per_sec=1000
changedate=datetime.date(2021, 5, 12) #date file changed FVwidth and decay time units

def find_files(subdir):
    if os.path.isdir(subdir):
        pattern = subdir+'*.pickle' #analyze set of files
        pickledir=subdir
    else:
        pattern=subdir+'.pickle' #analyze single file
        pickledir=os.path.dirname(pattern)+'\\'
    outfnames = sorted(glob.glob(pattern))
    if not len(outfnames):
        sys.exit('************ No files found *********')
    else:
        print('files=',outfnames)
        if '*' in pattern:
            pickledir=os.path.dirname(outfnames[0])+'\\'
    return outfnames,pickledir

def init_arrays(ps):
    ps.first_peak_end=int(ps.first_peak_end_fraction*ps.timepoints_per_trace)
    ps.latest_artifact=ps.artifact_window*ps.timepoints_per_trace 
    ps.artifactdecay=int(ps.artdecaytime/ps.dt)
    ps.FVwidth_pts=int(ps.FVwidth/ps.dt)
    ps.starttrace=0 #IO curve starts with the 1st trace
    ps.num_usable_traces=30 #IO curve consists of 1st 30 traces
    ps.num_usable_minutes=int(np.ceil(float(ps.num_usable_traces)/ps.tracesPerMinute)) #is this needed?
    ps.goodtraces=range(ps.starttrace,ps.starttrace+ps.num_usable_traces)
    print ("traces",ps.numtraces,"usable", ps.num_usable_traces,"start", ps.starttrace)
    #initialize arrays to hold measurements  
    ps.base=np.zeros((ps.num_usable_traces))
    ps.peak=np.zeros((ps.num_usable_traces))
    ps.peaktime=np.zeros((ps.num_usable_traces))
    ps.pospeak=np.zeros((ps.num_usable_traces))
    ps.pospeaktime=np.zeros((ps.num_usable_traces))
    ps.FVsize=np.zeros((ps.num_usable_traces))
    ps.amp=np.zeros((ps.num_usable_traces))
    ps.baseline_start=np.zeros((ps.num_usable_traces), dtype=int)
    ps.popspikestart=np.zeros((ps.num_usable_traces))        
    return

def analyze_IO(param_dict):
    params = Namespace(**param_dict)
    #
    #4. call PopSpikeAnalClass with the parameters, but do not generate plots
    pop_spike = psc.PopSpike(params) 
    pop_spike.read_datafile(question=False)
    init_arrays(pop_spike)
    pop_spike.findPopSpike()
    return pop_spike

def finish_IO(ps,start,end):
    ps.popspikenorm=np.zeros((ps.num_usable_minutes))
    ps.popspikeminutes=np.zeros((ps.num_usable_minutes))
    for sample in range(ps.num_usable_minutes):
        ps.popspikenorm[sample]=np.nanmean(ps.amp[sample*ps.tracesPerMinute:(sample+1)*ps.tracesPerMinute])	
        ps.popspikeminutes[sample]=(ps.tracetime[ps.starttrace+sample*ps.tracesPerMinute])/60	#60 sec per minute
    ps.current=np.linspace(start,end,len(ps.popspikenorm))

def showSummaryPlot(ps):
    from matplotlib import pyplot
    pyplot.ion()
    fig,axes=pyplot.subplots()
    fig.canvas.set_window_title('IO curve '+ps.params.exper)
    axes.plot(ps.current,ps.popspikenorm,'b.')
    axes.set_xlabel('Current (units)')
    axes.set_ylabel ('Mean Pop Spike (mV)')
    return    

def saveData(ps,new_outdir): 
    ps.params.exper=ps.params.exper+'_IO'
    ps.params.outputdir=new_outdir

    #tracedict includes all traces that you potentially want to analyze by groups in GroupAvg
    #modify argsdict to include all the other bits of information you want, and change dictionary keys to be meaningful for you.
    tracedict={'amp':ps.amp,'peaktime':ps.peaktime,'FVsize':ps.FVsize,'popspikenorm':ps.popspikenorm,'popspikeminutes':ps.popspikeminutes, 
                'current':ps.current}
    #This writes the file for GrpAvg.py
    #it includes all the experiment parameters
    datadict = dict(trace=tracedict, #datadict name wont be accessible outside this program but keys are
            parameters=ps.params, anal_params=ps.anal_params) #first is key second is value
    print(ps.params.outputdir,ps.params.exper)
    outpickle = ps.params.outputdir+ps.params.exper + '.pickle'
    with open(outpickle, 'wb') as f:
        pickle.dump(datadict, f) #saves the contents of datadict in file "f"...       
    return
    

if __name__=='__main__':
    #ARGS='C:\\Users\\kblackw1\\Documents\\Python_Scripts\\ephys_anal\\PopSpike2\\ C:\\Users\\kblackw1\\Documents\\Python_Scripts\\ephys_anal\\PopSpike2\\ 0.1,1.5'
    all_plots=False
    try:
        commandline = ARGS.split() # ARGS="pickle_dir new_outputdir "
        do_exit = False
    except NameError: 
        commandline = sys.argv[1:]
        do_exit = True    

    pickledir=commandline[0] #diirectory or single file
    if len(commandline)>2:
        new_outdir=commandline[2] #output directory
    else:
        new_outdir=commandline[0]
    start,end=[float(s) for s in commandline[1].split(',')] #start and end current
    pickle_files,pickledir=find_files(pickledir)
    lvm_dir=pickledir #specify directory of lvm here, if on a different computer

    for outfname in pickle_files:
        with open(outfname, 'rb') as f:
            params,_=get_params(f,outfname,do_exit) #make lvm_dir optional param?
            if not os.path.isdir(params['datadir']): #does previous datadir exist?
                params['datadir']=lvm_dir #try alternative datadir
            if params['exper'].endswith('_newb'):
               underscore=params['exper'].rfind('_')
               params['exper']=params['exper'][0:underscore]
            print('looking for',params['exper']+'.lvm','in',params['datadir'],', found:',glob.glob(params['datadir']+params['exper']+'.lvm'))
    
            pop_spike=analyze_IO(params)
            finish_IO(pop_spike,start,end)
            if all_plots:
                if len(pop_spike.problemtraces["tracenum"]):
                    pop_spike.plotProblemTraces()
                else:
                    print("!!!!!!No problem traces!!!!!!")
                pop_spike.plotGoodTraces()
            showSummaryPlot(pop_spike)
            if '*' in new_outdir:
                new_outdir=params['outputdir']
            saveData(pop_spike,new_outdir)
            
