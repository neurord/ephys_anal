#This program reads in all IF/IV data and creates SAS output.
#No possibility of using the parameter file to select or create subsets
#edit this to calculate IRneg and IRpos same as in MURI
#Also, calculate IFmax and IFhalfmax
import os
print os.getcwd()
#this is subdir for output data, relative to ~sarah
subdir="PickleIF/"
pattern = subdir+'*crv.pickle'
#no matter where I was when I started python, go to ~sarah to begin:
#home="/home/sarah/"
home="/home/avrama/EphysDataAnal/"
os.chdir(home)

import numpy as np
import sys  
from pprint import pprint as pp
import glob
import argparse
import pickle
import GrpPlotUtil as grp_utl
from matplotlib import pyplot
from matplotlib import gridspec
from scipy import optimize

printinfo=0
#It would be better to specify the  current injection values and then determine the points
#calcuate IRnegative as mean of values between points 0 and 7 (-450 pA to -200 pA)
IRneg_P1=0
IRneg_P2=5
#IRpos is mean between 40 pA and 100 pA
IRpos_P1=9
IRpos_P2=12
outfnames = sorted(glob.glob(pattern))
print "NUM FILES:", len(outfnames)

################# Select experiments that meet criteria
DATAS=[] #list of IF/IV data

for outfname in outfnames:
	with open(outfname) as fif:
		datadict = pickle.load(fif)
		DATAS.append(datadict)
		
if (len(DATAS)==0):
	print "no expers meet your criteria"
else:
        numexpers=len(DATAS)
        numcurrents=len(datadict['current'])
        #DATAS contain ALL time series data for ALL cells that are valid
        expername=[p['exper'] for p in  DATAS]
        width=[p['meanwidth'].x for p in  DATAS]
        height=[p['meanheight'].x for p in DATAS]
        rheobase=[p['rheobase'] for p in  DATAS]
        latency=[p['maxlat'] for p in  DATAS]
        meanahp=[p['mean_ahp'].x for p in  DATAS]
        spikes=[p['spike_count'] for p in  DATAS]
        currents=[p['current'] for p in  DATAS]
        spikefreq=[p['spike_freq'] for p in  DATAS]
        resistance =[p['resistance'] for p in  DATAS]
        IRneg=[np.mean(p['resistance'][IRneg_P1:IRneg_P2]) for p in  DATAS]
        IFmax=np.zeros(numexpers)
        IRpos=np.zeros(numexpers)
        IFhalfmax=np.zeros(numexpers)
        for exper in range(numexpers):
                if np.isnan(rheobase[exper]):
                        IRpos_end=IRpos_P2
                else:
                        IRpos_end=min(np.where(currents[exper]==rheobase[exper])[0][0],IRpos_P2)
                print exper,expername[exper],IRpos_end,rheobase[exper]
                IRpos[exper]=np.mean(resistance[exper][IRpos_P1:IRpos_end])
                IFmax[exper]=np.max(spikefreq[exper])
                #Find current at which spike _frequency_ is half of max
                if IFmax[exper]>0:
                        IFhalfmax[exper]=np.min(currents[exper][np.where(spikefreq[exper]>IFmax[exper]/2.0)])
                else:
                        IFhalfmax[exper]=np.nan
                #popt,pcov=optimize.curve_fit(grp_utl.sigmoid,currents[exper][np.where(spikes[exper]>0)],spikefreq[exper][np.where(spikes[exper]>0)])
                #print "other halfmax",IFhalfmax[exper],popt
                #IFhalfmax[exper]=popt[1]
        #
        ################ Write data for SAS
        SASoutput=np.column_stack((expername,rheobase,latency,width,height,meanahp,IRneg, IRpos,IFmax,IFhalfmax))
        SASheader="exper rheobase latency width height ahp IRneg IRpos IFmax IFhalfmax\n"
        f=open("IVIFforSAS.txt", 'w')
        f.write(SASheader)
        np.savetxt(f, SASoutput, fmt='%s', delimiter='   ')
        f.close()
os.chdir(home)
