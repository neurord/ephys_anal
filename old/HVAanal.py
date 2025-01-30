#This is the python program that reads in Calcium current traces from 2 Pulse protocol
#extracts baseline values, peak currents, steady state, and performs leak subtraction
#input is the full filename (including directory) of one of the traces, then python will find all the rest of the traces
#Outputs single file with each column containing a separate measure
#run the program from within python by typing: execfile('HVAanal.py')
#Possibly read in list of filenames, with start and end traces for epoch 1, 2 and 3, then do entire analysis in one go

import os
from igor import binarywave
from numpy import * # whats the asterisk doing again? importing All things in numpy? Yes.
from pylab import *
from matplotlib import pyplot
import sys  # if above is true, whats the dif here? ie why not from sys import *, or vice versa above?
from pprint import pprint as pp
import glob
from scipy import optimize

######################### Experiment specific parameters

#MUST BE SPECIFIED:
ExperName="082914_1"
ExperType=""
#minExperName is used to name output file, and igor waves containing output file
filenameending="_1_*_1_1p1.ibw"
IleakRamp5=0.0147046

#Should be OK from here
FileDir=ExperType+ExperName+"_Waves/"
pattern=FileDir+ExperType+"W"+ExperName+filenameending

print "Looking for files using: ", pattern

#some parameters common to ALL experiments, From RE analylsis
#Start fitting this amount of time after the prepulse begins
fitt1=0.0004
prepulsestart=0.05
#prepulsedur is actually end of 1 ms leak calc
prepulsedur=0.005
prepulseSSdur=0.001
#Measure baseline here to compare to pulse for leak calculation
pulseBaset1=0.045
#Next set of time points is to measure peak and ss of calcium currents
basestarttime=0.316
baseendtime=0.35
PeakStart=0.357
Peak2start=0.607
PeakDur=0.02
ssStart=0.544
ssEnd=ssStart+0.01
peakpoints=10
dt=50e-6       #50 microsecond sample rate
leakest=1
pulse2base=int(round(0.1/dt))
pulse2sst1=int(round(0.14/dt))
pulse2dur=int(round(0.004/dt))
pulse2=0
nA=1e9
ms=1e3
plotYN=1
#########################Below here, no parameters

def f(x,A,tau,B):
    return np.exp(-x/tau)*A+B

def f2(x,A1,tau1,A2,tau2,B):
    return np.exp(-x/tau1)*A1+np.exp(-x/tau2)*A2+B

def sortorder(fname):
    parts = fname.split('_')
    ans = int(parts[-3]), int(parts[-2])
    #print 'here:', fname, '->', ans
    return ans

filenames = glob.glob(pattern)
if (len(filenames)==0):
	print "Mistyped the filenames. Python found NO files:"
	pp(filenames)
else:

        filenames = sorted(filenames, key=sortorder)  #sorted is alphabetical unless integers given
                     #key is output of function "sortorder," using argument "filenames"
        #pp(filenames)

        #convert times into points
        prepulset1=prepulsestart+fitt1
        prepulset2=prepulsestart+prepulsedur
        prepulsept1=int(round(prepulset1/dt))
        prepulsept2=int(round(prepulset2/dt))
        prepulsepoints=int(round(prepulseSSdur/dt))
        pulseBaset1pt=int(round(pulseBaset1/dt))
        basestartpnt=int(round(basestarttime/dt))
        baseendpnt=int(round(baseendtime/dt))
        PeakStartpnt=int(round(PeakStart/dt))
        PeakEndpnt=int(round((PeakStart+PeakDur)/dt))
        Peak2startpnt=int(round(Peak2start/dt))
        Peak2endpnt=int(round((Peak2start+PeakDur)/dt))
        ssStartpnt=int(round(ssStart/dt))
        ssEndpnt=int(round(ssEnd/dt))
        fittime=linspace(fitt1,prepulsedur,(prepulsept2-prepulsept1+1))

        #initialize the arrays that will contain the data from all series and traces
        numtraces=len(filenames)
        Ihold=zeros(numtraces)
        Ileak5=zeros(numtraces)
        pulseSS=zeros(numtraces)
        IMin1=zeros(numtraces)
        IMin2=zeros(numtraces)
        Iss=zeros(numtraces)
        ss2peak=zeros(numtraces)
        mintime=zeros(numtraces)
        min2time=zeros(numtraces)
        peakRatio=zeros(numtraces)
        tracestring=[]
        base=zeros(numtraces)
        OptLeakEst=zeros(numtraces)
        tracedata=zeros((numtraces,20100))

        fig=pyplot.figure(figsize=(6,6))
        fig.canvas.set_window_title('Experiment '+ExperName)
        axes=fig.add_subplot(111)
    
        #loop through each trace from each series
    
        for fnum,eachfile in enumerate(filenames):
                data = binarywave.load(eachfile)      #read in data    
                trace=data['wave']['wData']
                tracedata[fnum]=trace
                #pp(data)             # optional - look at data array to see what else is there
                dtfile=data['wave']['wave_header']['hsA']  #verify that the dt we have is correct
                if (dt != dtfile):
                        print "Error, dt of file is different than expected"
                #
                parts = eachfile.split('_')
                tracestring.append(ExperName+'_'+parts[-3])
                #calculate various measures
                Ihold[fnum]=mean(trace[basestartpnt:baseendpnt])*nA
                IminPt=trace[PeakStartpnt:PeakEndpnt].argmin()+PeakStartpnt
                IMin1[fnum]=mean(trace[IminPt-peakpoints:IminPt+2*peakpoints])*nA
                mintime[fnum]=IminPt*dt
                IminPt2=trace[Peak2startpnt:Peak2endpnt].argmin()+Peak2startpnt
                IMin2[fnum]=mean(trace[IminPt2-peakpoints:IminPt2+2*peakpoints])*nA
                min2time[fnum]=IminPt2*dt
                Iss[fnum]=mean(trace[ssStartpnt:ssEndpnt])*nA
                #
                #Estimate the leak.  First, calc mean of last 1 ms, then curvefit the prepulse
                pulseSS[fnum]=mean(trace[prepulsept2-prepulsepoints:prepulsept2])*nA
                base[fnum]=mean(trace[pulseBaset1pt:pulseBaset1pt+prepulsepoints])*nA
                if leakest:
                    popt,pcov=optimize.curve_fit(f,fittime,trace[prepulsept1:prepulsept2+1])
                    Aopt,Tauopt,Bopt=popt
                    perr=np.sqrt(np.diag(pcov))
                    OptLeakEst[fnum]=Bopt*nA
                    Ileak5[fnum]=OptLeakEst[fnum]-base[fnum]
                else:
                    Ileak5[fnum]=pulseSS[fnum]-base[fnum]
                if pulse2:
                    pulse2ss=mean(trace[pulse2sst1:pulse2sst1+pulse2dur])
                    base2=mean(trace[pulse2base:pulse2base+pulse2dur])
                    Ileak10=(pulse2ss-base2)*nA
                #print out leak estimates to assess
                print tracestring[fnum], "leak ramp {0:.6f}, ss {1:.6f}, opt {2:.6f}".format(IleakRamp5,pulseSS[fnum]-base[fnum],OptLeakEst[fnum]-base[fnum]),
                if pulse2:
                    print "p2 {0:0.6f}".format(-0.5*Ileak10),
                if leakest:
                    print "Ierr {0:.6f}".format(perr[2]*nA), "tau {0:.4f} err {1:.5f}".format( Tauopt*ms,perr[1]*ms)
                #Optionally, plot the data to verify the results
                if plotYN:
                    endtime=len(trace)*dtfile
                    tracetime=arange(0,endtime,dt)
                    axes.plot(tracetime[0:len(trace)],trace*nA+fnum*0.1,label=parts[-3])
                    axes.plot(mintime[fnum],IMin1[fnum]+fnum*0.1,'ko') #label 1st min
                    axes.plot(min2time[fnum],IMin2[fnum]+fnum*0.1,'ro') #label 2nd min
                    axes.plot((ssStart+ssEnd)/2,Iss[fnum]+fnum*0.1,'bo') #label ss
                    axes.plot(pulseBaset1,base[fnum]+fnum*0.1, 'mo') #label baseline before prepulse
                    #label the leak subtraction values
                    if pulseSS[fnum]>base[fnum]:
                        axes.plot(prepulset2,pulseSS[fnum]+fnum*0.1,'go') #label the prepulse ss
                    else:
                        axes.plot(prepulset2,pulseSS[fnum]+fnum*0.1,'g*') #label the prepulse ss with *
                #
                #The leak estimates are not very good
                #replace with mean of first give ramp leak estimates (from RampAnal.py')
                Ileak50=IleakRamp5*10
                peakRatio[fnum]=(IMin2[fnum]-(Ihold[fnum]+Ileak50))/(IMin1[fnum]-(Ihold[fnum]+Ileak50))
                ss2peak[fnum]=(Iss[fnum]-(Ihold[fnum]+Ileak50))/(IMin1[fnum]-(Ihold[fnum]+Ileak50))
                print "Amp1 old {0:.6f} & new {1:.6f}".format(IMin1[fnum]-Ihold[fnum],IMin1[fnum]-(Ihold[fnum]+Ileak50)),
                print "ss2pk old {0:.6f} & new {1:.6f}".format((Iss[fnum]-Ihold[fnum])/(IMin1[fnum]-Ihold[fnum]), ss2peak[fnum])
                if pulse2:
                    axes.plot((prepulset2-prepulseSSdur),-0.5*Ileak10+fnum*0.1,'yo') 
                if leakest:
                    axes.plot(tracetime[prepulsept1:prepulsept2+1],f(fittime,Aopt,Tauopt,Bopt)*nA+fnum*0.1,'y')
        axes.legend(fontsize=8, loc='best')
        fig.canvas.draw()
        fig.show()

        #header contains wave names for igor
        header='trace'+ExperName+'  pkRatio'+ExperName+'  ss2pk'+ExperName+'  Ileak5'+ExperName+'  Ihold'+ExperName+'  Imin1'+ExperName+'   Imin2'+ExperName+'   Iss'+ExperName+'\n'

        outputdata=column_stack((tracestring,peakRatio,ss2peak,Ileak5,Ihold,IMin1,IMin2,Iss)) 
        outfname=ExperName+"HVAss2peak.txt"    
        #f=open(outfname,'w')
        #f.write(header)
        #savetxt(f, outputdata, fmt='%s', delimiter='     ')
        #f.close()
        outfname=ExperName+"traces.txt"
        outputdata=column_stack((tracetime,tracedata[0],tracedata[1],tracedata[2],tracedata[3],tracedata[4],tracedata[5],tracedata[6],tracedata[7],tracedata[8],tracedata[9],tracedata[10],tracedata[11],tracedata[12],tracedata[13],tracedata[14],tracedata[15],tracedata[16],tracedata[17]))
        f=open(outfname,'w')

        savetxt(f,outputdata)
        f.close()
