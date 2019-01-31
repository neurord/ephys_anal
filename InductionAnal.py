#Type something like this from within python: ARGS="2015_04_07_02"
#template: ARGS="filename"
#then type: execfile('InductionAnal.py')
#If not in python, type: python InductionAnal.py 031011_5Or
import os

#if running from windows, might need to execute the following command:
#os.getcwd()
#might need to use forward slashes in this command:
#pythondir="C:\Users\avrama\Documents\StochDiffcode\Python"
#os.chdir(pythondir)

import numpy as np
from matplotlib import pyplot
import sys
from pprint import pprint as pp
import glob
import pickle

#####parameters that you may want to tweak ###############
#If induction is not during the last pause in recording, you need to fix some of the code below
#Embedded in the head is the number of samples for one trace.  Best would be to read this values
sec_per_msec=0.001
plotYN=1

#Directory where data is located
datadir="E:\FIELDS\ValerieData//"
###################### below here, the only parameter is experiment naming convention #######################

try:
	fname = ARGS #in python: ARGS string
	do_exit = False
except NameError: #if you are not in python, read in filename 
	fname = sys.argv[1]
	do_exit = True

#Looks for file specified experiment name
filenamepattern=datadir+fname+".lvm"
print "Looking for files using: ", filenamepattern
filename = glob.glob(filenamepattern)
if (len(filename)==0):
	print "You mistyped the filenames. Python found no such files:"
	pp(filename)

#edit this to read from file
timepoints_per_trace=40000
#read in data file
data = np.loadtxt(filename[0], skiprows=24) 

#first column of datafile is time, second column is Vm
time=data[0:timepoints_per_trace,0]
tempVm=data[:,1]

dt=time[1]-time[0]

#reshape the data to put each trace in separate column
datalength=np.shape(tempVm)[0]
numtraces=datalength/timepoints_per_trace
Vm_traces=np.reshape(tempVm, (numtraces,timepoints_per_trace))

#plot the data
pyplot.ion()
fig,axes=pyplot.subplots()
axes.set_ylabel('Vm (mV)')
axes.set_xlabel('Time (sec)')
yoffset=0.1
xoffset=0
for tracenum in range(numtraces):
        axes.plot(time+tracenum*xoffset,Vm_traces[tracenum]+tracenum*yoffset)
fig,axes=pyplot.subplots()
axes.plot(data[:,0],data[:,1])
pyplot.show()
