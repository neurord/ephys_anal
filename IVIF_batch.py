#Python program to generate batch file for Ephys analysis processing

import os
from numpy import *
from pylab import *
from string import *
import glob
from pprint import pprint as pp
#######################################################
#Edit this part to indicate which files and which molecules
#Use glob on the .xml file
dir="IFcrv/"
suffix='ivif_Waves'
prninfo=1

python_prog="AnalyzeIV.py"
#######################################################

#1. Find all the files
pattern=dir+'*'+suffix
filenames = sorted(glob.glob(pattern))
pp(filenames)

#2. open file for writing
outfname='AnalyzeIV_batch.tmp'
f=open(outfname,'w')

#3. loop over all files and construct the processing line
for fullname in filenames:
    fname=fullname.split('/')[1]
    textline="python "+python_prog+" "+fname[0:find(fname,suffix)]+'\n'
    #output file for each input file
    f.write(textline) 
f.close()

        
