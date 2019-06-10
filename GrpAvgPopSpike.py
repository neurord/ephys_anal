#ARGS are optional, see ignore = for how they are used
#if using ARGS, from within python, type 
#ARGS="--sex M --age 21"
#execfile ('GrpAvgPopSpike.py')
#from outside python, type
#python GrpAvgPopSpike.py -s M -p 21 
import numpy as np
import sys  
from pprint import pprint as pp
import glob
import pickle
import os
print (os.getcwd())
from matplotlib import pyplot

#Change home (location of pickle files)
home="C:\Users\Sarah\Documents\Python Scripts\\"
os.chdir(home)

import pop_spike_utilities as psu
import GrpPlotUtil as grp_utl
#this is subdir for pickle files.  Uncomment the appropriate one:
#subdir="AlexData\\"
subdir="Pickle\\"
plot_individuals=0  #to plot PopSpike vs time for each experiment in a group
plot_corr=0 #to plot  correlation between LTP (at summarytime) and age or baseline epsp
print_info=0
colors=pyplot.get_cmap('viridis')
color_incr=32

#VARIABLE MAY CHANGE; 15min baseline plus 30min follow-up is 90sweeps:
summarytime=[25,30]
minimum_sweeps=45      # <<<<<<<<<<<<<<<<<< units are minutes
slope_std_factor=2     #<<<<<<<<<<<<<<<<<<< if slope +/- slope_std_factor*std includes 0, trace is valid, else slope too large
                       #Sarah actually allows a negative slope with potentiation and vice versa
slope_threshold=0.01  #0.008-0.01 give good no-stim control
nan_threshold=5     #do not use experiment if too many missing (nan) popspikes 
########## Specify separation variables for generating means to plot in Igor	
###For binary separation, list only one value; for > 2 values, must list each one
#if you list only one separation variable, you only get 2 windows in the graph
sepvarlist=[['region',['DM']],['theta',[0, 5.0, 10.5]],['drug',['none','nitr', 'dmso', 'sch23390']]]#,['age',[50]]],#,['sex',['F','Fe','M']]]#
#sepvarlist=[['theta',[0, 5.0, 10.5]],['drug', ['none']]]#,['age',[50]]]
pattern = subdir+'*.pickle'
outfnames = sorted(glob.glob(pattern))
print ("NUM FILES (before selection):", len(outfnames))

try:
	commandline = ARGS.split() #in python: define space-separated ARGS string
	do_exit = False
except NameError: #if you are not in python, read in filename and other parameters
	commandline = sys.argv[1:]
	do_exit = True

specify_params=psu.parse_args(commandline,do_exit,1)

################# Select experiments that meet criteria
DATAS = [] # list will contain all trace data from experiments which are not ignored
PARAMS = [] #list will contain all the parameters, so we can create SAS output
ANAL=[]
BAD=[]
fig,axes=pyplot.subplots(2,1)
fig.canvas.set_window_title('Problems')
axes[0].set_ylabel('popspike-bad slope')
axes[1].set_ylabel('popspike-nan')
axes[1].set_xlabel('Time (min)')
nan_color=0
nancount=[]

for outfname in outfnames:
    with open(outfname) as f:
        datadict = pickle.load(f)
    if print_info:
        print ("file read:", datadict['parameters'],"baseline slope", round(datadict['trace']['slope'],6))
    exper_param = datadict['parameters']
    numnan=sum(np.isnan(datadict['trace']['popspikenorm']))
    ignore1 = ((specify_params.sex and specify_params.sex != exper_param.sex) or 
                  (specify_params.age is not None and specify_params.age >= exper_param.age) or
                  (specify_params.maxage is not None and specify_params.maxage <= exper_param.age) or
                  (specify_params.drug and specify_params.drug != exper_param.drug) or
                  (specify_params.region and specify_params.region != exper_param.region) or
                  (specify_params.theta and specify_params.theta != exper_param.theta))
    ignore2 = ((len(datadict['trace']['popspikenorm'])<minimum_sweeps) or
                  (np.abs(datadict['trace']['slope'])>slope_threshold) or
                  #-slope_std_factor*datadict['trace']['slope_std']
                  numnan>nan_threshold) #do not use trace if too many missing popspikes                                           

    if ignore2 and not ignore1:
        if len(datadict['trace']['popspikenorm'])<minimum_sweeps:
            print ("!!!NOT ENOUGH goodtraces", datadict['parameters'].exper, len(datadict['trace']['popspikenorm']))
        if np.abs(datadict['trace']['slope'])>slope_threshold:
            print ("!!!BAD baseline", datadict['parameters'].exper,  "slope u,s", round(datadict['trace']['slope'],6),round(datadict['trace']['slope_std'],6))
            axes[0].plot(datadict['trace']['popspikeminutes'],datadict['trace']['popspikenorm'],'.', label=datadict['parameters'].exper)
            BAD.append(datadict)
        #text=raw_input('continue? (y/n)')
    if ignore1 or ignore2:
        if print_info:
            print ("***** ignoring:", datadict['parameters'])
        next
    else:
        if np.isnan(datadict['trace']['popspikenorm']).any():
            print ("########## np.nan detected", exper_param.exper, datadict['anal_params'])
            axes[1].plot(datadict['trace']['popspikeminutes'],datadict['trace']['popspikenorm'],'+', label=datadict['parameters'].exper)
            nan_index=np.argwhere(np.isnan(datadict['trace']['popspikenorm']))
            axes[1].plot(datadict['trace']['popspikeminutes'][nan_index],np.ones((len(nan_index))),'o', label=datadict['parameters'].exper)
        if datadict['trace']['popspikenorm'][-1]==0.0:
            print ("@@@@@@@@@ check popSpikeAnal for", exper_param.exper, datadict['anal_params'])
        #print 'OK: {}'.format(exper_param)
        DATAS.append(datadict['trace'])
        PARAMS.append(datadict['parameters'])
        ANAL.append(datadict['anal_params'])
        #
        if datadict['trace']['popspikeminutes'][1]-datadict['trace']['popspikeminutes'][0]<0.8:
            print(datadict['parameters'].exper,datadict['anal_params']['artdecay'], datadict['anal_params']['FVwidth'],datadict['trace']['popspikeminutes'])
axes[0].legend(fontsize=8, loc='best')
axes[1].legend(fontsize=8, loc='best')
fig.canvas.draw()
if len(BAD):
    print ("&&&&&&&&&&&&& Bad Baselines")
    from scipy import optimize
    badexper=[p['parameters'].exper for p in BAD]
    badps=[p['trace']['popspikenorm'] for p in BAD]
    badx=[p['trace']['popspikeminutes'] for p in BAD]
    badslope=[p['trace']['slope'] for p in BAD]
    for i in range(len(badexper)):
        baseline_end=np.min(np.where(badx[i]>0))
        validbasepts=~np.isnan(badps[i][0:baseline_end])
        popt,pcov=optimize.curve_fit(psu.line,badx[i][validbasepts],badps[i][validbasepts])
        Aopt,Bopt=popt
        Astd,Bstd=np.sqrt(np.diag(pcov))
        validbasepts=~np.isnan(badps[i][5:baseline_end])
        popt,pcov=optimize.curve_fit(psu.line,badx[i][validbasepts],badps[i][validbasepts])
        Aopt,shortBopt=popt
        Astd,shortBstd=np.sqrt(np.diag(pcov))
        print ("bad baseline", badexper[i], "slope",round(Bopt,5), "+/-", round(Bstd,5), "10min slope", round(shortBopt,6), round(shortBstd,5))
#Calculate average over all data that meets criteria
if (len(DATAS)==0):
	print ("no expers meet your criteria")
else:
        #DATAS and datas contain ALL time series data for ALL cells that are valid
        ##can also add some summary measures for each cell, e.g. mean baselime, mean normpeakamp at a few time points
        Sex=[p.sex for p in PARAMS]
        Age=[p.age for p in PARAMS]
        region=[p.region for p in PARAMS]
        Drug=[p.drug for p in PARAMS]
        theta=[p.theta for p in PARAMS]
        exper=[p.exper for p in PARAMS]
        slope=[p['slope'] for p in DATAS]
        slope_std=[p['slope_std'] for p in DATAS]
        ps_means=[p['PS_mean'] for p in DATAS]
        fv_means=[p['FV_means'] for p in DATAS]
        artifacthresh=[p['artifactthresh'] for p in ANAL]
        FVwidth=[p['FVwidth'] for p in ANAL]
        artdecay=[p['artdecay'] for p in ANAL]
        noisethresh=[p['noisethresh'] for p in ANAL]
        dict_grp=DATAS
        ################ Write all of valid data for SAS: parameters and 5 min means of popspikenorm (ps_mean) and FVnorm (fv_means)
        SASoutput=np.column_stack((exper,Sex,Age,Drug,region,theta,ps_means, fv_means))
        SASheader="exper Sex Age drug region theta normpopspike norm_fv\n"
        f=open("PARAMSforSAS.txt", 'w')
        f.write(SASheader)
        np.savetxt(f, SASoutput, fmt='%s', delimiter='   ')
        f.close()
        Experheader="exper Sex Age drug region theta artifacthresh FVwidth artdecay noisethresh\n"
        Exper_list=np.column_stack((exper,Sex,Age,Drug,region,theta,artifacthresh,FVwidth, artdecay, noisethresh))
        f=open("Experimentlist.txt", 'w')
        f.write(Experheader)
        np.savetxt(f, Exper_list, fmt='%s', delimiter='   ')
        f.close()
        #
        ######Separate out data into multiple arrays based on categorical variables
        for sepnum in range(len(sepvarlist)):
            sepvar=sepvarlist[sepnum][0]
            sepvalue=sepvarlist[sepnum][1]
            if sepnum==0:
                #print "******first separation:",sepvar,', value=', sepvalue, 
                data_to_separate=DATAS #holds popspikenorm
                params_to_separate=PARAMS
                dict_grp,param_grp,=grp_utl.separate(data_to_separate,sepvar,sepvalue,params_to_separate)
                #print ', total exp:',len(DATAS),'; grp',sepvar,"=",sepvalue[0],'has',len(dict_grp[0]),'; grp', sepvar,"NE",sepvalue[0],'has',len(dict_grp[1])
            else:
                #print "  ***additional separations:", sepvar, ', value=', sepvalue
                data_to_separate=dict_grp
                dict_grp=[]
                params_to_separate=param_grp
                param_grp=[]
                for num in range(np.shape(data_to_separate)[0]):
                    #print "    before sep of", sepvar, ", group:", num, "#exp=", np.shape(data_to_separate[num])[0]
                    if len(data_to_separate[num])>1:
                        tempdata,tempparam=grp_utl.separate(data_to_separate[num],sepvar,sepvalue,params_to_separate[num])
                        for i in range(len(tempdata)):
                                dict_grp.append(tempdata[i])
                                param_grp.append(tempparam[i])
                                if len(sepvalue)==1:
                                        sep_phrase=" NE "+str(sepvalue[0])
                                else:
                                        sep_phrase=" = "+str(sepvalue[i])
                                #print "       after ",'split; grp',sepvar,sep_phrase,'has', len(tempdata[i])
                    elif len(data_to_separate[num])==1:
                            dict_grp.append(data_to_separate[num])
                            param_grp.append(params_to_separate[num])
                            #print "       one exper in group - appending, no further splitting"
                    #else:
                            #print "       zero exper in group - no further splitting"
        #
        numgroups=len(dict_grp)
        avgpopspikenorm_grp=[]
        stderrpopspikenorm_grp=[]
        newcount_grp=[]
        minutes_grp=[]
        filenm=[]
        paramgrpnum=[]
        pnum=0
        for groupnum in range(numgroups):
            #print "\ncalculate average for: groupnum", groupnum,", numexpers=", len(dict_grp[groupnum]),
            if len(dict_grp[groupnum])>1:
                tempavg,tempstd,tmpcount= grp_utl.exp_avg(dict_grp[groupnum],'popspikenorm')
                avgpopspikenorm_grp.append(tempavg)
                stderrpopspikenorm_grp.append(tempstd/np.sqrt(len(dict_grp[groupnum])))
                newcount_grp.append(tmpcount)
                tempavg,tempstd,tmpcount= grp_utl.exp_avg(dict_grp[groupnum],'popspikeminutes')
                minutes_grp.append(tempavg)
                filenm.append(grp_utl.construct_filename(sepvarlist,param_grp[groupnum]))
                paramgrpnum.append(pnum)
                pnum=pnum+1
                #print ', filename=',filenm[-1]
            else:
                if len(dict_grp[groupnum]):
                    avgpopspikenorm_grp.append(dict_grp[groupnum][0]['popspikenorm'])
                    stderrpopspikenorm_grp.append(np.zeros(len(dict_grp[groupnum][0]['popspikenorm'])))
                    newcount_grp.append(np.ones(len(dict_grp[groupnum][0]['popspikenorm'])))
                    minutes_grp.append(dict_grp[groupnum][0]['popspikeminutes'])
                    filenm.append(grp_utl.construct_filename(sepvarlist,param_grp[groupnum]))
                    paramgrpnum.append(pnum)
                    pnum=pnum+1
                    print (', filename=',filenm[-1])
                else:
                    pnum=pnum+1
                    numgroups=numgroups-1
                    #print ', no output file for this group'
        ####### Write output
        print ("&&&&&&&&&&&&& Summary Data")
        for gnum in range(numgroups):
            #construct output filename that tells you which experiments, e.g. Valid PSP
            f=open(filenm[gnum]+".txt",'w')                      
            outputdata=np.column_stack((minutes_grp[gnum], newcount_grp[gnum],
                                        100*avgpopspikenorm_grp[gnum], 100*stderrpopspikenorm_grp[gnum]))
            header="time "+filenm[gnum]+"count "+filenm[gnum]+"normpsAVG "+filenm[gnum]+"normpsSE "
            f.write(header+"\n")
            np.savetxt(f,outputdata, fmt='%7.5f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
            f.close()
            #print some summary information, 25 min = 40 samples
            firstpt=np.min(np.where(minutes_grp[gnum]>summarytime[0]))
            lastpt=np.max(np.where(minutes_grp[gnum]<summarytime[1]))
            print ('Grp',filenm[gnum], 'at', summarytime,'minutes: ', np.mean(avgpopspikenorm_grp[gnum][firstpt:lastpt]),' n=',newcount_grp[gnum][0])
            print ('SD', np.mean(stderrpopspikenorm_grp[gnum][firstpt:lastpt])*np.sqrt(newcount_grp[gnum][0]))
        ########## End of data separation
        ######## Plot results
        grp_utl.plot_groups(avgpopspikenorm_grp,stderrpopspikenorm_grp,minutes_grp,newcount_grp,filenm,sepvarlist,param_grp,paramgrpnum)

print (str(len(DATAS))+'experiments met criteria.')

if plot_individuals:
  for gnum in range(numgroups):
    grp_utl.plot_onegroup(dict_grp[paramgrpnum[gnum]],param_grp[paramgrpnum[gnum]], filenm[gnum])       

if plot_corr:
    min25_to_30=1
    min55_to_60=2
    sampletime=min25_to_30
    grp_utl.plot_corr(dict_grp,param_grp,paramgrpnum,numgroups,filenm,sampletime)

