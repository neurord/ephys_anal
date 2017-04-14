#This program reads in EPSP traces from plasticity experiments
#extracts baseline values, peak PSP, etc
#input is data and exp# for a single experiment, for example "12-Jun-2014_SLH006"
#input arguments by defining ARGS as a single string (dividing arguments with spaces)
#outputs include one file with each column containing a separate measure (for IGOR)
#outputs include dictionary housing exp features used to select like exps for analysis (temp, age, etc)
#Can use Zbyszek's python program to read in IF or IV curves and extract spike measures
#run the program from within python by typing: execfile('PSPanalSA.py')
#ARGS="20-Aug-2015_SLH001 M 29 heat nodrug 5.4 12 A2a+ MSN non 0 soma APs" for example

import os
import pylab
from igor import binarywave
import numpy as np
from matplotlib import pyplot
import sys
from pprint import pprint as pp
import glob
import argparse
import pickle

#####parameters I may want to tweak on viewing plot of data###############

PSPstart=0.262#generally good here-------- 
PSPend=0.285 #&&&&<<<---------MANUALLY define ... generally 0.265-0.28
peak2exists=True #&&&&<<<----MANUALLY define, True if compound/multi EPSP


pp2exists= False#ALWAYS FALSE... throwback to when test pulses were 50ms-spaced pairs
PP2start=0.3      #------------------------------------MANUALLY define

basestarttime=0.15  
baseendtime=0.25
hyperstart=0.045    #10 - 60 ms, last 20 ms is steady state          
hyperend=0.065
dt=100e-6       
plotstart=0 # change to "baseendtime/dt" to zoom in on EPSP in plots
Iaccess=50e-12 # -50 pA current injected to monitor series R

##########################################################################
######################### Experiment type specific parameters

#Descriptive variables used in subsequent program to select similar data to average:
parser = argparse.ArgumentParser()
parser.add_argument('experiment', type=str, help = 'give exp name, for example: 12-Jun-2014_SLH006')
parser.add_argument('--no-graphs', '-g', dest='graphs', default=True, action='store_false') # -g optional and not position-defined (its absence OK too)
parser.add_argument("sex", type=str, choices=["M","F"],help="male M or female F")
parser.add_argument("age", type=int, help="animal age in days")
parser.add_argument("bathtemp", type=str, choices=["heat","RT"],help="heat or RT")
parser.add_argument("drug", type=str, choices=["nodrug","NorBNI","Nomi","CGP","other"],help="nodrug,NorBNI,Nomi,CGP,other")
parser.add_argument("Rtip", type=float,help="pipette tip resistance")
parser.add_argument("indtime", type=int,help="minutes between break-in and TBS starting")
parser.add_argument("genotype", type=str, choices=["D1+", "A2a+", "D1-", "A2a-", "wt"],help="D1+ A2a+ (or-) for cre lines, or wt")
parser.add_argument("celltype", type=str, choices=["MSN", "FSI", "other"],help="MSN, FSI, or other cell type")
parser.add_argument("lightresp", type=str, choices=["LR", "non"],help="was cell light responsive? (LR or non)")
parser.add_argument("lightlevel", type=int, help="what light level was used? 0-100")
parser.add_argument("depol", type=str, choices=["soma", "soma0", "opto"],help="depolarization during TBS: soma, soma0, or opto")
parser.add_argument("TBSAP", type=str, choices=["APs","noAPs"],help="during TBS, APs or noAPs")
keydict={"sex":["M","F"],"bathtemp":["heat","RT"],"drug":["nodrug","NorBNI","Nomi","CGP","Naloxone","other"],"genotype":["D1+", "D2+", "D1-", "D2-", "wt"],"celltype":["MSN", "FSI", "other"],"lightresp":["LR", "non"],"depol":["soma", "soma0", "opto"],"TBSAP":["APs","noAPs"]}

with open('choicedict.txt', "wb") as f:# opens the file, for writing 
	pickle.dump(keydict, f) # writes the keydict to an external file named choicedict 
try:
	commandline = ARGS.split() #in python: define space-separated ARGS string
	do_exit = False
except NameError: #undefined variable (in this case ARGS)
	commandline = sys.argv[1:]
	do_exit = True

try:
	args = parser.parse_args(commandline) # maps arguments (commandline) to choices, and checks for validity of choices. see args.sex,age,etc.. printed below
#if arguements are mapped incorrectly, python wants to exit, but the next line says "don't", instead check whether we are in python (do_exit=False) then don't exit, just give us a warning
except SystemExit:
	if do_exit:
		raise # raise the exception above (SystemExit) b/c none specified here
	else:
		raise ValueError('invalid ARGS')

#print experiment characteristics to double check unconstrained entries
print "sex={}".format(args.sex)
print "age={}".format(args.age)
print "bathtemp={}".format(args.bathtemp)
print "drug={}".format(args.drug)
print "Rtip={}".format(args.Rtip)
print "indtime={}".format(args.indtime)
print "genotype={}".format(args.genotype)
print "celltype={}".format(args.celltype)
print "lightresp={}".format(args.lightresp)
print "lightlevel={}".format(args.lightlevel)
print "depol={}".format(args.depol)
print "TBSAP={}".format(args.TBSAP)

#Looks for files using pattern to find all files w/in a single exper
outputDir="Pickle/"
filenameending="_SLH_1_*_*_1p1.ibw"
inputDir="PatchData/"
#To re-analyze slope from all experiments, 1st use glob on FileDir
#2nd, read in pickel file for the parameters
#3rd, add in analysis of slope
#4th, resave with slope
#Can the computer detect existance of 2nd peak?
#Implement stuff avrama had tried, and plot which has 2nd peak and which doesn't, for sarah to view
FileDir=inputDir+args.experiment+"_Waves/"
parts=args.experiment.split('-')
pattern=FileDir+"W"+parts[0]+"_"+parts[1]+"_"+parts[2].split("_")[0]+filenameending
print "Looking for files using: ", pattern
filenames = glob.glob(pattern)
if (len(filenames)==0):
	print "You mistyped the filenames. Python found no such files:"
	pp(filenames)


#Below puts your files from common exper in order by PGF and trace
def sortorder(fname):
    parts = fname.split('_')
    ans = int(parts[-3]), int(parts[-2])
    #print 'here:', fname, '->', ans
    return ans
filenames = sorted(filenames, key=sortorder)#sorted is alphabetical unless integers given
                     #key is output of function "sortorder," using argument "filenames"
#pp(filenames)

#convert times into points
basestartpnt=basestarttime/dt
baseendpnt=baseendtime/dt
PSPstartpnt=PSPstart/dt
PSPendpnt=PSPend/dt
PP2startpnt=PP2start/dt          #<----------------------------------------------
hyperstartpnt=hyperstart/dt
hyperendpnt=hyperend/dt

#additional measures: tracetime
#initialize the arrays that will contain the data from all series and traces

RMP=[]

peakvm=[]
peaktime=[]
peakamp=[]

peak2=[]
peak2time=[]
peak2amp=[]

pp2=[]
pp2time=[]
pp2amp=[]

hypervm=[]
hyperdelta=[]

AccR=[]

tracestring=[]



fig=pyplot.figure(figsize=(6,6))
fig.canvas.set_window_title('Experiment '+args.experiment)
axes=fig.add_subplot(111)
    
#loop through each trace from each series
tempRMP=[]
tempAccess=[]
goodtraces=0
badcount=0
reallybad=0
baselineVm=0
for i,filename in enumerate(filenames):
        data = binarywave.load(filename)   
        trace=data['wave']['wData']
        Vm = np.mean(trace[basestartpnt:baseendpnt]) 
        tempRMP.append(Vm)
	access=(tempRMP[i]-np.mean(trace[hyperstartpnt:hyperendpnt]))/Iaccess/1e6 #converts access units to megaohms
	tempAccess.append(access)

        if (i==9): 
            baselineVm=np.mean(tempRMP) #this is mean over traces
	    baselineAccess=np.mean(tempAccess)

	if i>=10:
	    if (abs(Vm-baselineVm)/abs(baselineVm))> 0.2 or (abs(access-baselineAccess)/abs(baselineAccess))>.4: #20% baseline change or 40% access change
                reallybad = True
            else:
                reallybad = False

        bad  = Vm > -0.060 # <------- -0.060 
        if bad or reallybad:
            badcount += 1
        elif badcount <= 10: #< 10.  varying would make +- lax criteria -----------
            badcount = 0
            goodtraces = i+1	
	#print "VM", Vm, "Badcount", badcount,'# goodtraces is ', goodtraces, ', meaning ', (goodtraces/2)-5, ' minutes of follow-up.'

print '# goodtraces is ', goodtraces, ', meaning ', (goodtraces/2)-5, ' minutes of follow-up.'

for fileindex,filename in enumerate(filenames[:goodtraces]):
        data = binarywave.load(filename)      #read in data  
        trace=data['wave']['wData']#wave is NOT x data, but wData is y data.. type "data" into python to explore binary raw data!
        # pp(data) 
	#wave is key housing dict with key wave_header housing dict with key hsA... which is key linked to dt
        dtfile=data['wave']['wave_header']['hsA']  #verify that the dt we have is correct
        if dt != dtfile:
            raise ValueError("Error, dt of file is different than expected")
        
#calculate various measures.
        RMP.append(np.mean(trace[basestartpnt:baseendpnt]))

	maxvm_i = trace[PSPstartpnt:PSPendpnt].argmax() + PSPstartpnt #argmax gives index (point number) of maximum, within slice
        maxvm = trace[maxvm_i]
	peakvm.append(maxvm)
        peakamp.append(peakvm[fileindex]-RMP[fileindex])
	maxvm2_i = trace[PSPendpnt:].argmax() + PSPendpnt
	peaktime.append(maxvm_i*dt) #<------maxvm_i gives point number of maximum.. why multiply by dt??

	if pp2exists == True:
	    maxvmpp2_i = trace[PP2startpnt:].argmax() + PP2startpnt
            maxvmpp2 = trace[maxvmpp2_i]
	    pp2.append(maxvmpp2)
            pp2amp.append(pp2[fileindex]-RMP[fileindex])
	    pp2time.append(maxvmpp2_i*dt) 

        if peak2exists == True:
            maxvm2 = trace[maxvm2_i]
	    peak2.append(maxvm2)
            peak2time.append(maxvm2_i*dt)
	    peak2amp.append(peak2[fileindex]-RMP[fileindex])


	hypervm.append(np.mean(trace[hyperstartpnt:hyperendpnt]))
        hyperdelta.append(RMP[fileindex]-hypervm[fileindex])
	AccR.append(hyperdelta[fileindex]/Iaccess)
    	parts = filename.split('_')
        tracestring.append(parts[-3]+"_"+parts[-2]+"   ") #extra spaces forces column_stack to use more significant figures when converting floats to string. There should be a better way<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ?

        #print filename.split('/')[1],tracestring[fileindex],RMP[fileindex],peakvm[fileindex]
        #Optionally, plot the data to verify the results
        endtime=len(trace)*dtfile
        tracetime=np.arange(0,endtime,dt) #array of points between 0 and endtime stepping by dt
        if args.graphs:     
            x=PSPstart
            #print x
            y_array=np.arange(-.085,-.06)
            x_array=x*np.ones(y_array.shape)
            if fileindex<10:
	    	axes.plot(tracetime[plotstart:],trace[plotstart:],label=fileindex, color=(0,fileindex*.1, fileindex*.1))
                axes.plot(tracetime[PSPstartpnt],trace[PSPstartpnt],'rs',ms=10)
		axes.plot(tracetime[PSPendpnt],trace[PSPendpnt],'rd',ms=10)
		axes.plot(tracetime[PP2startpnt],trace[PP2startpnt],'gd',ms=10)
            elif fileindex>(goodtraces-11):
		axes.plot(tracetime[PSPstartpnt],trace[PSPstartpnt],'rs',ms=10)
		axes.plot(tracetime[PSPendpnt],trace[PSPendpnt],'rd',ms=10)
	        axes.plot(tracetime[plotstart:],trace[plotstart:],label=fileindex, color=((fileindex+11-goodtraces)*.1,0,(fileindex+11-goodtraces)*.1))
            	axes.plot(tracetime[PP2startpnt],trace[PP2startpnt],'gd',ms=10)
        #if fileindex==2:
        #    break    
        #axes.plot(peaktime,peakvm[fileindex]-RMP[fileindex]+fileindex*0.001,'ko') #label the peak
	#axes.plot(peak2time,peak2[fileindex]-RMP[fileindex]+fileindex*0.001,'ko') #label the peak

axes.legend(fontsize=8, loc='best')
fig.canvas.draw()
#comment out the next line to suppress graphs
fig.show()

#calculating change from baseline
deltavm=np.mean(RMP[0:10])-RMP # in mV, subtraction from baseline avg RMP 
normpeakamp = peakamp/np.mean(peakamp[0:10])
normpeaktime=peaktime/np.mean(peaktime[0:10])
if peak2exists == True:
    normpeak2amp=peak2amp/np.mean(peak2amp[0:10]) 
    normpeak2time=peak2time/np.mean(peak2time[0:10])
else:
    normpeak2amp=[]
    normpeak2time=[]
normAccR=AccR/np.mean(AccR[0:10])
normRMP=RMP/np.mean(RMP[0:10])

index=np.arange(len(peakamp))
if args.graphs:
	#before saving data formatted for Igor, just evaluate the results of this experiment
	fig=pyplot.figure(figsize=(10,10))
	fig.canvas.set_window_title('Summary '+args.experiment)
	

	axes=fig.add_subplot(411)# says of 4 plots, we're dealing with 1st column 2nd row
	axes.plot(index,normpeakamp,'rp',label='peak1, % change')
	if peak2exists == True:
	    axes.plot(index,normpeak2amp,'.',label='peak2, % change')
	#axes.axis([0,len(index),0,2])
	axes.legend(fontsize=10, loc='best')

	axes=fig.add_subplot(412)# 4 plots can happen, we're dealing with 1st column top row
	axes.plot(index,peakamp,'rp',label='peak1 (V)')
	if peak2exists == True:
	    axes.plot(index,peak2amp,'.',label='peak2 (V)') #if error, may need to place this in  an "if exists" context
	#axes.axis([0,len(index),0,.03])
	axes.legend(fontsize=10, loc='best')

	axes=fig.add_subplot(413)# says of 4 plots, we're dealing with 1st column 3rd row
	axes.plot(index,normAccR-1,'-',label='access-R change vs baseline')
	#axes.plot(index,normpeaktime-1,'@',label='% change in time to 1st peak')
	#axes.plot(index,deltavm*10,'*',label='deltaVm*10') # *10 so I can see it on graph... 
	#axes.axis([0,len(index),-0.2,.2])
	axes.legend(fontsize=10, loc='best')

	axes=fig.add_subplot(414)# 4th of 4 plots in one column
	axes.plot(index,RMP,'*',label='RMP (mV)')
	#axes.axis([0,len(index),-0.055,-.095])
	axes.legend(fontsize=10, loc='best')

fig.canvas.draw()
#comment out the next line to suppress graphs
pyplot.show()

#modify tracedict to include all traces that you want, modify argsdict to include all the other bits of information you want, and change dictionary keys to be meaningful for you.
tracedict={'trace':tracestring,'RMP':RMP,'normRMP':normRMP,'deltavm':deltavm,'peakamp':peakamp,'pp2amp':pp2amp,'peak2amp':peak2amp,
'normpeakamp':normpeakamp,'normpeak2amp':normpeak2amp,'peaktime':peaktime,'peak2time':peak2time,'normpeaktime':normpeaktime,
'normpeak2time':normpeak2time,'AccR':AccR, 'normAccR':normAccR, 'goodtraces':goodtraces}
#one dictionary for experiment type which has only 2 values
#argsdict={'temp':args.bathtemp,'age':args.age,'sex':args.sex}
#another dictoinary for experiment type which has a bunch of values, and you want to use > or <

datadict = dict(trace=tracedict, #datadict name wont be accessible outside this program but keys are
		parameters=args) #first is key second is value

#This writes the file for GrpAvg.py
#it includes all the experiment parameters (args)
outfname=outputDir+args.experiment+"wcDATA"
outpickle = outfname + '.pickle'
with open(outpickle, 'w') as f:
	pickle.dump(datadict, f) #saves the contents of datadict in file "f"... name datadict meaningless outside but contents (parameters, etc) retain meaning.

#The next lines creates a text file for reading into igor (not needed to read into GrpAvg.py)
#, but doesn't include experiment type information
#header contains wave names for igor   #add goodminutes variable here perhaps? Where? Want to access it in next program
#Need to clean this up. <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if peak2exists == True:
	header=args.experiment+'index '+args.experiment+'trace '+args.experiment+'base '+args.experiment+'peak '+args.experiment+'amp '+args.experiment+'hyper'+args.experiment+'amp2\n'#+args.experiment+'PPamp2\n'
	outputdata=np.column_stack((index,tracestring,RMP,peakvm,peakamp,hyperdelta,peak2amp))#,pp2amp)) 
else:
	header=args.experiment+'index '+args.experiment+'trace '+args.experiment+'base '+args.experiment+'peak '+args.experiment+'amp '+args.experiment+'hyper\n'#+args.experiment+'PPamp2\n'
	outputdata=np.column_stack((index,tracestring,RMP,peakvm,peakamp,hyperdelta))#,pp2amp)) 

  
#<<<<<<<<<<<<<<<<< WHAT DOES THIS CHUNK DO? It opens and writes the file for "typical PSP time course (not trace)" 
with open(outfname+".txt",'w') as f: # dont use "close" with "with" b/c its built in... "with" just makes sure the file finish/closes after being opened. The variable f is set equal to the product of the open() clause
	f.write(header) 
	np.savetxt(f, outputdata, fmt='%s', delimiter='     ')



