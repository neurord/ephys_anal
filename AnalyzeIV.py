#To run from within python: 1. ARGS="14-Mar-2015_SLH002" 2. execfile('AnalyzeIV.py')
# start and end of current duration, and time to wait until steady state, should be specified here
#need to fix spike height, currently it is the overshoot
from loader import Measurement
import numpy as np
import pickle
import sys
#Do this whenever you change some function definitions
#import loader; reload(loader); from loader import Measurement

#The directory names where input data are stored (subdir) and pickle files are written (outdir)
experType="crv"
dirnameEnding="ivif_Waves"
subdir='IFcrv/'
outdir='PickleIF/'
#bad_traces=([5,1],[5,2],[5,3],[5,4],[5,5])
#bad_traces=([5,6])
#duration of current injection
injectdur=0.4

try:
    commandline = ARGS.split()  #from within python, type ARGS="14-Mar-2015_SLH002", optionally include number of IVseries
    do_exit = False
except NameError: 
    commandline=sys.argv[1:]   #from outside of python; e.g. "10-Dec-2014_SLH009"
    do_exit = True

Exper=commandline[0]
if len(commandline)>1:
    IVnum=int(commandline[1])
else:
    IVnum=4
dirname=subdir+Exper+dirnameEnding
#default: assumes IVseries is 4 and IF is anything else
#default IF current begins at 200 pA and goes up by 20
#default IV current begins at -500 pA and goes up by 50
waves1 = Measurement(dirname,IF=(40e-12, 20e-12),IV=(-450e-12, 50e-12),time=0.8,IVseries=IVnum,)
#waves1 = Measurement(dirname,IF=(40e-12, 20e-12),IV=(-450e-12, 50e-12),time=0.8,IVseries=IVnum,bad_sweep=bad_traces)

#### This starts our code, assuming measurements1 exists and works
spike_count=waves1.spike_count
latency=waves1.spike_latency
current=waves1.injection
width=waves1.spike_width
meanwidth=waves1.mean_spike_width
########### Why no access to spike_height (only means)
#height=waves1.spike_height
meanheight=waves1.mean_spike_height
if np.size(np.where(spike_count>0)):
    rheobase=current[np.min(np.where(spike_count>0))]
    maxlat=np.max(latency[np.isfinite(latency)])
else:
    rheobase=np.nan
    maxlat=np.nan
ahp=waves1.spike_ahp
mean_ahp=waves1.mean_spike_ahp
deltaV=waves1.response
IR=np.zeros(len(deltaV))
spike_freq=np.zeros(len(spike_count))

series_flag=1
for i in range(len(deltaV)):
    if spike_count[i]==0:
        IR[i]=deltaV[i][0]/current[i]/1e6         #Input resistance in Units of Megaohms
        spike_freq[i]=0
    else:
        spike_freq[i]=spike_count[i]/(injectdur-latency[i])
        IR[i]=np.nan
    if IR[i]<0:
        print "whoops, wrong currents, IR negative for", Exper, "IVnum=",IVnum 
        series_flag=0       

if series_flag==1:
    print "writing pickle file"
    #plot IF curves with spike peak, AHP min marked
    IFdict={'exper':Exper, 'current': current, 'DeltaV': deltaV, 'spike_count':spike_count,'latency':latency,'resistance':IR,
            'width': width,'ahp':ahp,'meanheight':meanheight,'meanwidth':meanwidth,'rheobase':rheobase,'maxlat':maxlat,'mean_ahp':mean_ahp,'spike_freq':spike_freq}

    #Need to specify or construct outfname
    outfname=Exper+experType
    outpickle = outdir+outfname + '.pickle'
    #with open(outpickle, 'w') as f:
    #    pickle.dump(IFdict, f) #saves the contents of datadict in file "f"... name datadict meaningless outside but 
