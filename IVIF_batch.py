#Python program to generate batch file for Ephys analysis processing

from numpy import *
import glob
#######################################################
#Edit this part to indicate which files 
#Use glob on the .xml file
dir="IFcrv/"
suffix='ivif_Waves'

python_prog="AnalyzeIV.py"
#######################################################

#1. Find all the files
pattern=dir+'*'+suffix
filenames = sorted(glob.glob(pattern))

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

        
