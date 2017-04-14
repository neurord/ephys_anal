#This is the python program that reads in test pulse traces
#extracts baseline values, peak PSP, etc
#input is the full filename (including directory) of one of the traces, then python will find all the rest of the traces
#Outputs single file with each column containing a separate measure
#Can use Zbyszek's python program to read in IF or IV curves and extract spike measures
#run the program from within python by typing: execfile('PSPanalSA.py')

import os
from igor import binarywave
from numpy import * # whats the asterisk doing again? importing All things in numpy? Yes.
from pylab import *
from matplotlib import pyplot
import sys  # if above is true, whats the dif here? ie why not from sys import *, or vice versa above?
from pprint import pprint as pp
import glob

######################### Experiment specific parameters

#MUST BE SPECIFIED:
ExperName="092114_3"
ExperType="ramp"
#minExperName is used to name output file, and igor waves containing output file
filenameending="_1_*.ibw"

#Should be OK from here
FileDir=ExperType+ExperName+"_Waves/"
pattern=FileDir+ExperType+ExperName+filenameending

print "Looking for files using: ", pattern

#some parameters common to ALL experiments, From RE analylsis
#perhaps the base needs to be changed - go from the max to account for leak?
basestarttime=0.028  
baseEndtime=0.046
RampStart=0.05   
RampVslope=0.5  #mV/ms = V/sec
RampEnd=0.31
PeakStart=0.100
PeakEnd=0.300
maxEnd=0.15
dt=50e-6       #50 microsecond sample rate
#estimate leak from the ramp around -60 mV:
#V=-65 mV at t=0.08 sec and -55 mv at t=0.1 sec)
rampleak1=int(round(0.08/dt))
rampleak2=int(round(0.10/dt))
DeltaV=-10e-3      #10 mV
rampPts=5
#########################Below here, no parameters

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
        basestartpnt=basestarttime/dt
        baseEndpnt=baseEndtime/dt
        PeakStartpnt=PeakStart/dt
        PeakEndpnt=PeakEnd/dt
        maxEndpt=maxEnd/dt

        #initialize the arrays that will contain the data from all series and traces
        numtraces=len(filenames)
        Ihold=zeros(numtraces)
        IMin=zeros(numtraces)
        Amp=zeros(numtraces)
        IMax=zeros(numtraces)
        Amp2=zeros(numtraces)
        tracestring=[]
        #-0.5e9 converts to nA and leak due to 5 mV (not 10 mV) Vm change
        RampLeak=zeros(numtraces)

        fig=pyplot.figure(figsize=(6,6))
        fig.canvas.set_window_title('Experiment '+ExperName)
        axes=fig.add_subplot(111)
    
        #loop through each trace from each series
    
        for fileindex,eachfile in enumerate(filenames):
                parts = eachfile.split('_')
                tracestring.append(ExperName+'_'+parts[-3]+"_"+parts[-2]+"      ") #extra spaces forces column_stack to use more significant figures when converting floats to string.  There should be a better way
                data = binarywave.load(eachfile)      #read in data    
                trace=data['wave']['wData']
                LeakVsTime=zeros(len(trace))
                #pp(data)             # optional - look at data array to see what else is there
                dtfile=data['wave']['wave_header']['hsA']  #verify that the dt we have is correct
                if (dt != dtfile):
                        print "Error, dt of file is different than expected"
                #calculate various measures
                Ihold[fileindex]=mean(trace[basestartpnt:baseEndpnt])
                IMin[fileindex]=min(trace[PeakStartpnt:PeakEndpnt])
                mintimePt=trace[PeakStartpnt:PeakEndpnt].argmin()+PeakStartpnt
                mintime=mintimePt*dt
                maxtimePt=trace[PeakStartpnt:maxEndpt].argmax()+PeakStartpnt
                maxtime=maxtimePt*dt
                IMax[fileindex]=mean(trace[maxtimePt-10:maxtimePt+10])
                Amp[fileindex]=Ihold[fileindex]-IMin[fileindex]
                #Calculate slope of I vs V from linear part:
                DeltaI=mean(trace[rampleak1-rampPts:rampleak1+rampPts])-mean(trace[rampleak2-rampPts:rampleak2+rampPts])
                LeakSlope=DeltaI/DeltaV
                #-0.5e-9 converts leak into value in response to 5 mV depol
                RampLeak[fileindex]=-0.5e9*(DeltaI)
                #Calculate Ramp current versus time
                endtime=len(trace)*dtfile
                tracetime=arange(0,endtime-dt,dt)
                LeakVsTime=LeakSlope*RampVslope*(tracetime-RampStart)+Ihold[fileindex]
                LeakVsTime[0:int(RampStart/dt)]=Ihold[fileindex]
                LeakVsTime[int(RampEnd/dt):]=Ihold[fileindex]
                Amp2[fileindex]=LeakVsTime[mintimePt]-IMin[fileindex]
                LeakSubtrRamp=trace-LeakVsTime
                print tracestring[fileindex], "ramp Leak {0:.6f}".format( RampLeak[fileindex] )
                #Optionally, plot the data to verify the results
                #axes.plot(tracetime,trace+fileindex*0.1e-9,label=fileindex)
                axes.plot(tracetime,LeakSubtrRamp+fileindex*0.1e-9,label=fileindex)
                #axes.plot(tracetime,LeakVsTime)
                axes.plot(0.12,LeakVsTime[mintimePt]+fileindex*0.1e-9,'ko') #label the min
                axes.plot(mintime,IMin[fileindex]-LeakVsTime[mintimePt]+fileindex*0.1e-9,'b*') #label the min
        axes.legend(fontsize=8, loc='best')
        fig.canvas.draw()
        fig.show()

        endLeak=min(5,len(RampLeak))
        print "endleak", tracestring[endLeak-1],"Mean Leak", np.mean(RampLeak[0:endLeak])

        #header contains wave names for igor
        header='trace'+ExperName+'    Amp2_'+ExperName+'    Amp'+ExperName+'    Leak'+ExperName+'\n'

        outputdata=column_stack((tracestring,Amp2*1e9,Amp*1e9,RampLeak)) 
        outfname=ExperName+"RampAmp.txt"    
        f=open(outfname,'w')
        f.write(header)
        savetxt(f, outputdata, fmt='%s', delimiter='     ')
        f.close()
        #outputdata=column_stack((tracetime,trace,LeakVsTime,LeakSubtrRamp))
        #outfname=ExperName+"RampLeakSubtr.txt"
        #f=open(outfname,'w')
        #savetxt(f,outputdata)
        #f.close()
