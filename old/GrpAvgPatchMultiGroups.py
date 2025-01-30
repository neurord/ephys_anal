#ARGS are optional, see ignore = for how they are used
#if using ARGS, from within python, type 
#ARGS="-bt heat -Rt 4.5 -Rtmax 6.8 -geno D1+ -cell MSN"
#execfile ('GrpAvgPatchMultiGroups.py')
#from outside python, type
#python GrpAvg.py -s M -p 21 -bt RT

#Needs to be updated to use plot_groups from grp_utl
#Needs to be updated (along with PSPanalSA.py) to use common argparser

import os
print os.getcwd()
#home is location of this file and GrpPlotUtil
home="/home/sarah/"
#home="/home/avrama/EphysDataAnal/"
os.chdir(home)

#this is subdir for output data (from PSPanalSA.py), relative to home
subdir="Pickle/"
pattern = subdir+'*wcDATA.pickle_new'

import numpy as np
import sys  
from pprint import pprint as pp
import glob
import argparse
import pickle
import GrpPlotUtil as grp_utl
from matplotlib import pyplot
from matplotlib import gridspec

printinfo=0
#VARIABLE MAY CHANGE:
meanstart=40   #trace number corresponding to 15 min after induction
meanend=50     #trace number corresponding to 20 min after induction
minsweeps=50        # minimum sweeps for completed experiment
normpeakthresh=4	#if higher than this, then AP must have occured, exclude from mean
########## Specify separation variables for generating means to plot in Igor	
###### This first param in the list will be used to plot two different groups, averaged over the remaining params
llval=0
sepvarlist=[['lightresp',['LR']],['drug',['nodrug','NorBNI','Nomi','CGP']],['lightlevel',[llval]]]##,['compound',0]]#['genotype',['D1']],['cre',['+']],

parser = argparse.ArgumentParser()
parser.add_argument('--no-graphs', '-g', dest='graphs', default=True, action='store_false') # just type -g
parser.add_argument("--sex", '-s', type=str, choices=["M","F"],help="male M or female F")
parser.add_argument("--age", '-p', type=int, help="minimum animal age in days")
parser.add_argument("--bathtemp", '-bt',type=str, choices=["heat","RT"],help="heat or RT")
parser.add_argument("--drug", '-d', type=str, choices=["nodrug","NorBNI","Naloxone","other"],help="nodrug,NorBNI,Naloxone,other")
parser.add_argument("--Rtip", '-Rt', type=float,help="minimum pipette tip resistance")
parser.add_argument("--Rtipmax", '-Rtmax', type=float,help="max pipette tip resistance")
parser.add_argument("--indtime", '-indtime',type=float,help="max acceptable minutes between break-in and TBS")
parser.add_argument("--indtimemin", '-indmin', type=float,help="min acceptable minutes to wait on stable baseline")
parser.add_argument("--genotype", '-geno',type=str, choices=["D1+", "D1-", "wt","A2a+","A2a-"],help="D1+ A2a+ (or-) for cre lines, or wt")
parser.add_argument("--celltype", '-cell',type=str, choices=["MSN", "FSI", "other"],help="MSN, FSI, or other cell type")
parser.add_argument("--lightresp", '-li',type=str, choices=["LR", "non"],help="was cell light responsive? (LR or non)")
parser.add_argument("--lightlevel", '-levmin',type=float, help="minimum light level, 0-100")
parser.add_argument("--lightlevelmax", '-levmax',type=float, help="maximum light level, 0-100")
parser.add_argument("--depol", '-dep',type=str, choices=["soma", "soma0", "opto"],help="depolarization during TBS: soma, soma0, or opto")
parser.add_argument("--TBSAP", '-APs',type=str, choices=["APs","noAPs"],help="during TBS, APs or noAPs")
actual = parser.parse_args()

outfnames = sorted(glob.glob(pattern))
print "NUM FILES:", len(outfnames)

try:
	commandline = ARGS.split(" ") #within python: define variable ARGS as a string housing space-separated arguments
	print "ARGS =", ARGS, "commandline=", commandline
 	do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
	commandline = sys.argv[1:]
	print "commandline =", commandline
	do_exit = True

try:
	args = parser.parse_args(commandline) # maps arguments (commandline) to choices, and checks that they match
			
except SystemExit:
	if do_exit:
		raise
	else:
		raise ValueError('invalid ARGS')

################# Select experiments that meet criteria
DATAS = [] # list will contain all trace data from experiments which are not ignored
PARAMS = [] #list will contain all the parameters, so we can create SAS output

for outfname in outfnames:
	with open(outfname) as f:
		datadict = pickle.load(f)
        #print datadict['parameters'].experiment , datadict['traces]['goodtraces']
	#pp(datadict)
	#pp(datadict['parameters'])
	actual = datadict['parameters']
	#print 'actual', actual
	ignore = ((args.sex and args.sex != actual.sex) or #I set an option, and its != whats in datadict
                  (args.age and args.age >= actual.age) or
                  (args.bathtemp and args.bathtemp != actual.bathtemp)or
		  (args.drug and args.drug != actual.drug) or
                  (args.Rtip and args.Rtip > actual.Rtip) or
                  (args.Rtipmax and args.Rtipmax < actual.Rtip) or
                  (args.lightlevel is not None and args.lightlevel >= actual.lightlevel) or
                  (args.lightlevelmax is not None and args.lightlevelmax <= actual.lightlevel) or
		  (args.indtimemin and args.indtimemin >= actual.indtime) or
		  (args.indtime and args.indtime <= actual.indtime) or
                  (args.genotype and args.genotype != actual.genotype) or
                  (args.celltype and args.celltype != actual.celltype) or
                  (args.lightresp and args.lightresp != actual.lightresp) or
                  (args.depol and args.depol != actual.depol) or
                  (args.TBSAP and args.TBSAP != actual.TBSAP) or
	          len(datadict['trace']['RMP'])<minsweeps)
	if ignore:
		next 
		#print "ignoring:",datadict['parameters'].experiment ,"   of length", len(datadict['trace']['RMP']#['traces']['goodtraces']
                #print "ignoring:", datadict['experiment'], "goodtraces", datadict['trace']['goodtraces']
                #text=raw_input('continue? (y/n)')
                #print "ignoring:",datadict['parameters'].experiment, " of len", len(datadict['trace']['RMP']), actual
	else:
		#print 'OK: {}'.format(actual) 
		DATAS.append(datadict['trace'])
		PARAMS.append(datadict['parameters'])
				
if (len(DATAS)==0):
	print "no expers meet your criteria"
else:
        #DATAS and datas contain ALL time series data for ALL cells that are valid
        ##can also add some summary measures for each cell, e.g. mean baselime, mean normpeakamp at a few time points
        Sex=[p.sex for p in PARAMS]
        Age=[p.age for p in PARAMS]
        Batht=[p.bathtemp for p in PARAMS]
        Drug=[p.drug for p in PARAMS]
        Rtip=[p.Rtip for p in PARAMS]
        Indt=[p.indtime for p in PARAMS]
        Geno=[p.genotype for p in PARAMS]
        Cellt=[p.celltype for p in PARAMS]
        Lightr=[p.lightresp for p in PARAMS]
        Lightlev=[p.lightlevel for p in PARAMS]
        Depol=[p.depol for p in PARAMS]
        TBSAP=[p.TBSAP for p in PARAMS]
        exper=[p.experiment for p in PARAMS]
        compound=[p.compound for p in PARAMS]
        dict_grp=DATAS

        ################ Write all of valid data for SAS: parameters and 5 min means of normpeak, normnewpeak, normslope
        timeframes=10
        samples_per_frame=10
        normpeakSAS=np.zeros((len(PARAMS),timeframes))
        normnewpeakSAS=np.zeros((len(PARAMS),timeframes))
        normslopeSAS=np.zeros((len(PARAMS),timeframes))
        basePSP=np.zeros(len(PARAMS))
        baseslope=np.zeros(len(PARAMS))
        for i in np.arange(len(DATAS)):	
                DATAStemp=DATAS[i]['normpeakamp']
                DATAStemp[DATAStemp > normpeakthresh] = np.nan
                for j in np.arange(timeframes):
                        DATAScount=(DATAStemp[j*samples_per_frame:(j+1)*samples_per_frame])		
                        normpeakSAS[i,j]=np.nansum(DATAScount)/np.sum(-np.isnan(DATAScount))#count if NOT(-) not a number
                DATAStemp=DATAS[i]['normslope']
                for j in np.arange(timeframes):
                        DATAScount=(DATAStemp[j*samples_per_frame:(j+1)*samples_per_frame])		
                        normslopeSAS[i,j]=np.nansum(DATAScount)/np.sum(-np.isnan(DATAScount))#count if NOT(-) not a number
                DATAStemp=DATAS[i]['normnewpeakamp']
                for j in np.arange(timeframes):
                        DATAScount=(DATAStemp[j*samples_per_frame:(j+1)*samples_per_frame])		
                        normnewpeakSAS[i,j]=np.nansum(DATAScount)/np.sum(-np.isnan(DATAScount))#count if NOT(-) not a number
                DATAStemp=DATAS[i]['peakamp']
                basePSP[i]=np.mean(DATAStemp[0:10])
                DATAStemp=DATAS[i]['slope']
                baseslope[i]=np.nanmean(DATAStemp[0:10])
       #
        SASoutput=np.column_stack((exper,Sex,Age,Batht,Drug,Rtip,Indt,Geno,Cellt,Lightr,Lightlev,Depol,TBSAP,basePSP,compound, normpeakSAS,normnewpeakSAS,baseslope, normslopeSAS))
        SASheader="exper Sex Age Batht Drug Rtip Indt Geno Cellt Lightr Lightlev Depol TBSAP basePSP compound normpeak normnewpeak baseslope normslope\n"
        f=open("PARAMSforSAS.txt", 'w')
        f.write(SASheader)
        np.savetxt(f, SASoutput, fmt='%s', delimiter='   ')
        f.close()

        ######Separate out data into multiple arrays based on categorical variables
        #print "beginning separation"
        for sepnum in range(len(sepvarlist)):
            sepvar=sepvarlist[sepnum][0]
            sepvalue=sepvarlist[sepnum][1]
            if sepnum==0:
		if printinfo:
                	print "******first separation:",sepvar,', value=', sepvalue, 
                data_to_separate=DATAS #holds 
                params_to_separate=PARAMS
                dict_grp,param_grp,=grp_utl.separate(data_to_separate,sepvar,sepvalue,params_to_separate)
		if printinfo:
                	print ', total exp:',len(DATAS),'; grp',sepvar,"=",sepvalue[0],'has',len(dict_grp[0]),'; grp', sepvar,"NE",sepvalue[0],'has',len(dict_grp[1])
            else:
		if printinfo:
                	print "  ***additional separations:", sepvar, ', value=', sepvalue
                data_to_separate=dict_grp
                dict_grp=[]
                params_to_separate=param_grp
                param_grp=[]
                for num in range(np.shape(data_to_separate)[0]):
		    if printinfo:
                    	print "    before sep of", sepvar, ", group:", num, "#exp=", np.shape(data_to_separate[num])[0]
                    if len(data_to_separate[num])>1:
                        tempdata,tempparam=grp_utl.separate(data_to_separate[num],sepvar,sepvalue,params_to_separate[num])
                        for i in range(len(tempdata)):
                                dict_grp.append(tempdata[i])
                                param_grp.append(tempparam[i])
                                if len(sepvalue)==1:
                                        sep_phrase=" NE "+str(sepvalue[0])
                                else:
                                        sep_phrase=" = "+str(sepvalue[i])
				if printinfo:
                                	print "       after ",'split; grp',sepvar,sep_phrase,'has', len(tempdata[i])
                    elif len(data_to_separate[num])==1:
                            dict_grp.append(data_to_separate[num])
                            param_grp.append(params_to_separate[num])
			    if printinfo:
                            	print "       one exper in group - appending, no further splitting"
                    else:
			    if printinfo:
                            	print "       zero exper in group - no further splitting"
        #
        #AVERAGE over the experiments in each group
        #
        numgroups=len(dict_grp)
        avgnormpeakamp_grp=[]
        stderrnormpeakamp_grp=[]
        avgnormslope_grp=[]
        stderrnormslope_grp=[]
        avgnormnewpeak_grp=[]
        stderrnormnewpeak_grp=[]
        newcount_grp=[]
        filenm=[]
        paramgrpnum=[]
        minutes_grp=[]
        pnum=0
        for groupnum in range(numgroups):
		if printinfo:
                	print "\ncalculate average for: groupnum", groupnum,", numexpers=", len(dict_grp[groupnum]),
                if len(dict_grp[groupnum])>1:
			if printinfo:
                        	print "avg > 1: groupnum", groupnum,"numexpers", len(data_grp[groupnum])
                        tempavg,tempstd,tmpcount= grp_utl.exp_avg(dict_grp[groupnum],'normpeakamp',normpeakthresh)
                        avgnormpeakamp_grp.append(tempavg)
                        stderrnormpeakamp_grp.append(tempstd/np.sqrt(len(dict_grp[groupnum])))
                        tempavg,tempstd,tmpcount= grp_utl.exp_avg(dict_grp[groupnum],'normnewpeakamp',normpeakthresh)
                        avgnormnewpeak_grp.append(tempavg)
                        stderrnormnewpeak_grp.append(tempstd/np.sqrt(len(dict_grp[groupnum])))
                        newcount_grp.append(tmpcount)
                        tempavg,tempstd,tmpcount= grp_utl.exp_avg(dict_grp[groupnum],'normslope',normpeakthresh)
                        avgnormslope_grp.append(tempavg)
                        stderrnormslope_grp.append(tempstd/np.sqrt(len(dict_grp[groupnum])))
                        #construct output filename that tells you which experiments, e.g. Valid PSP
                       	filenm.append(grp_utl.opto_filename(sepvarlist,param_grp[groupnum],llval))
                        paramgrpnum.append(pnum)
                        pnum=pnum+1
			if printinfo:
                        	print ', filename=',filenm[-1]
                else:
                        if len(dict_grp[groupnum]):
				if printinfo:
                         		print "avg == 0 or 1: groupnum", groupnum,"numexpers", len(data_grp[groupnum])
                                avgnormpeakamp_grp.append(dict_grp[groupnum][0]['normpeakamp'])
                                avgnormnewpeak_grp.append(dict_grp[groupnum][0]['normnewpeakamp'])
                                avgnormslope_grp.append(dict_grp[groupnum][0]['normslope'])
                                newcount_grp.append(np.ones(len(dict_grp[groupnum][0]['normpeakamp'])))
                                stderrnormpeakamp_grp.append(np.zeros(len(dict_grp[groupnum][0]['normpeakamp'])))
                                stderrnormnewpeak_grp.append(np.zeros(len(dict_grp[groupnum][0]['normnewpeakamp'])))
                                stderrnormslope_grp.append(np.zeros(len(dict_grp[groupnum][0]['normslope'])))
                                filenm.append(grp_utl.opto_filename(sepvarlist,param_grp[groupnum],llval))
                                paramgrpnum.append(pnum)
                                pnum=pnum+1
				if printinfo:
                                	print ', filename=',filenm[-1]
                        else:
                                pnum=pnum+1
                                numgroups=numgroups-1
				if printinfo:
                                	print ', no output file for this group'
        ######### Prepare output
	tot=0
        for gnum in range(numgroups):
		if printinfo:
                	print "prep output: groupnum", gnum,"numexpers", len(dict_grp[paramgrpnum[gnum]]), 
                if len(dict_grp[paramgrpnum[gnum]])>0:
                        index=np.arange(len(avgnormpeakamp_grp[gnum]))
                        minutes_grp.append(index/2.0)
			if printinfo:
                        	print filenm[gnum]
                        f=open(filenm[gnum]+".txt",'w')                 
                        outputdata=np.column_stack((minutes_grp[gnum], newcount_grp[gnum],
                                                    avgnormpeakamp_grp[gnum], stderrnormpeakamp_grp[gnum],
                                                    avgnormnewpeak_grp[gnum], stderrnormnewpeak_grp[gnum],
                                                    avgnormslope_grp[gnum], stderrnormslope_grp[gnum])) 
                        header="time "+filenm[gnum]+"count "+filenm[gnum]+"normpeakAVG "+filenm[gnum]+"normpeakSE "+filenm[gnum]+"normnewpeakAVG "+filenm[gnum]+"normnewpeakSE "+filenm[gnum]+"normslopeAVG "+filenm[gnum]+"normslopeSE "
                        f.write(header+"\n")
                        np.savetxt(f,outputdata, fmt='%7.4f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
                        f.close()
			tot=tot+newcount_grp[gnum][0]
                        print '       Grp',filenm[gnum], '15-20 minutes out: ', np.mean(avgnormpeakamp_grp[gnum][meanstart:meanend]),' n=',newcount_grp[gnum][0],'tot=',tot
                else:
                        print "no outputfile"
                #
                ########## End of data separation, averaging, and output
        ######## Begin plotting data.  Traces are separated on the basis of first two separation variables
        pyplot.ion()
	fig=pyplot.figure(figsize=(15,15))
	fig.canvas.set_window_title(str(len(DATAS))+' selected Exps, averaged')
        numwindows=4
        gs = gridspec.GridSpec(numwindows, 2) #n rows, 2columns
        maxdur=0
        for num in range(len(avgnormpeakamp_grp)):
                if minutes_grp[num][-1]>maxdur:
                        maxdur=minutes_grp[num][-1]
        graphnum=[]
        for num in range(len(avgnormpeakamp_grp)):
                if (getattr(param_grp[paramgrpnum[num]][0],sepvarlist[1][0])==sepvarlist[1][1][0]) and (getattr(param_grp[paramgrpnum[num]][0],sepvarlist[0][0])==sepvarlist[0][1][0]):
                        graphnum.append(0)
                if (getattr(param_grp[paramgrpnum[num]][0],sepvarlist[1][0])!=sepvarlist[1][1][0]) and (getattr(param_grp[paramgrpnum[num]][0],sepvarlist[0][0])==sepvarlist[0][1][0]):
                        graphnum.append(1)
                if (getattr(param_grp[paramgrpnum[num]][0],sepvarlist[1][0])==sepvarlist[1][1][0]) and (getattr(param_grp[paramgrpnum[num]][0],sepvarlist[0][0])!=sepvarlist[0][1][0]):
                        graphnum.append(2)
                if (getattr(param_grp[paramgrpnum[num]][0],sepvarlist[1][0])!=sepvarlist[1][1][0]) and (getattr(param_grp[paramgrpnum[num]][0],sepvarlist[0][0])!=sepvarlist[0][1][0]):
                        graphnum.append(3)
        #
        for gridnum in range(numwindows):
                axesa=fig.add_subplot(gs[gridnum,0])
                axesb=fig.add_subplot(gs[gridnum,1])
                axesa.axhline(1)
                axesa.axis([0,maxdur,.4,2.0])
                axesa.set_ylabel('normpeakamp')
                axesb.axhline(1)
                axesb.axis([0,maxdur,.4,2.0])
                axesb.set_ylabel('normnewpeakamp')
                if gridnum==0:
                        axesa.set_xlabel('Time (min)')
                        axesb.set_xlabel('Time (min)')
                for num in range(len(avgnormpeakamp_grp)):
                        if graphnum[num]==gridnum:
                                axesa.errorbar(minutes_grp[num],avgnormpeakamp_grp[num],stderrnormpeakamp_grp[num],label=filenm[num]+',n='+str(newcount_grp[num][0]))
                                axesb.errorbar(minutes_grp[num],avgnormslope_grp[num],stderrnormslope_grp[num],label=filenm[num]+',n='+str(newcount_grp[num][0]))
                axesa.legend(fontsize=10, loc='best')
        fig.canvas.draw()
        pyplot.show()
        #execfile('SarahGraphs.py')
print str(len(DATAS))+'experiments met criteria.'
os.chdir(home)
