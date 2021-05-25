import numpy as np
from matplotlib import pyplot
from matplotlib import gridspec
from scipy import optimize
import pop_spike_utilities as psu

######## This function calculates averages across experiments
def exp_avg(datas, datatype,somethreshold=999e6):
    if datatype=='none':
        alldatavalues=datas
        #print "exp_avg: datatype=none, single array passed"
    else:
        alldatavalues = [celldata[datatype] for celldata in datas]
    ln = max(len(celldatavalues) for celldatavalues in alldatavalues)
    total = np.zeros(ln)
    count = np.zeros(ln)
    sumsquares = np.zeros(ln)
  
    for celldatavalues in alldatavalues:
        for i,datavalue in enumerate(celldatavalues):
            if datavalue < somethreshold and ~np.isnan(datavalue):
                total[i]+=datavalue
                count[i]+=1
                sumsquares[i]+=datavalue**2
			#else:
				#print (datatype,datavalue)
    mn=total/count #cannot use np.nanmean because the lengths of each array in alldatavalues are not identical
    meansquares=sumsquares/count
    stdev=meansquares-mn**2
    return mn,stdev,count

def separate(in_data,sepvar,sepval,in_params):
    #make sep_data and sep_params have as many sub arrays as length of sepval?
    separate_data=[[],[]]
    sep_params=[[],[]]
    for i in range(len(sepval)-2): 
        sep_params.append([])
        separate_data.append([])
    #print "sep_params",sep_params, "sepval", sepval, sepval[0]
    for i in range(len(in_data)): 
        success=0
        if sepvar=='genotype':
            attribute=getattr(in_params[i],sepvar)[0:-1]
        else:
            if sepvar=='cre':
                attribute=getattr(in_params[i],'genotype')[-1]
            else:
                attribute=getattr(in_params[i],sepvar)
        if len(sepval)==1:
            if isinstance(sepval[0], str):
               	#two classes, either equal or not equal to specified param
               	if (str(attribute) == sepval[0]):
                   	separate_data[0].append(in_data[i])
                   	sep_params[0].append(in_params[i])
                   	success=1
               	else:
                   	separate_data[1].append(in_data[i])
                   	sep_params[1].append(in_params[i])
                   	success=1
            else:
            	if (attribute > sepval[0]):
					#print "sepval", sepval[0], '<' ,attribute
                    separate_data[0].append(in_data[i])
                    sep_params[0].append(in_params[i])
                    success=1
            	else:
                    #rint "sepval", sepval[0], '>=', attribute
                     separate_data[1].append(in_data[i])
                     sep_params[1].append(in_params[i])
                     success=1				
        else: 
            #multiple classes - each value has to be specified
            for j in range(len(sepval)):
                if (attribute == sepval[j]):
                    sep_params[j].append(in_params[i])
                    separate_data[j].append(in_data[i])
                    success=1
            #print "separate on", sepvar, attribute
        if success==0:
            print ("!!!!!!!!!!data",i,"  not assigned to group:",in_params[i])
        #else:
                #print "data",i," assigned", in_params[i]
    return separate_data,sep_params

def plot_groups(avg_grp,stderr_grp,minutes_grp,count,filenm,sepvarlist,paramgrp,paramgrpnum):
    window_criterion=2 #make this 1 to get 4 windows with 2 separation variables
    #maxdur=minutes_grp[0][-1]
    pyplot.ion()
    fig=pyplot.figure(figsize=(10,15))
    fig.canvas.set_window_title('Group Average of Normalized PopSpike')
    if len(sepvarlist)>window_criterion:
        numwindows=4
    else:
        numwindows=2
    gs=gridspec.GridSpec(numwindows,1)
    #Assign graph number based on first two separation variables
    graphnum=[]
    for num in range(len(avg_grp)):
        if len(sepvarlist)<=window_criterion:
            if (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[0][0])==sepvarlist[0][1][0]):
                graphnum.append(0)
            else:
                graphnum.append(1)
        if len(sepvarlist)>window_criterion:
            if (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[1][0])==sepvarlist[1][1][0]) and (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[0][0])==sepvarlist[0][1][0]):
                graphnum.append(0)
            if (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[1][0])==sepvarlist[1][1][0]) and (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[0][0])!=sepvarlist[0][1][0]):
                graphnum.append(1)
            if (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[1][0])!=sepvarlist[1][1][0]) and (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[0][0])==sepvarlist[0][1][0]):
                graphnum.append(2)
            if (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[1][0])!=sepvarlist[1][1][0]) and (getattr(paramgrp[paramgrpnum[num]][0],sepvarlist[0][0])!=sepvarlist[0][1][0]):
                graphnum.append(3)
    #plot the data
    for gridnum in range(numwindows):
        axes=fig.add_subplot(gs[gridnum,0])
        axes.axis([-15,60,.6,1.6])
        axes.axhline(1)
        axes.set_ylabel('normPopSpike')
        axes.set_xlabel('Time (min)')
        for num in range(len(avg_grp)):
            if graphnum[num]==gridnum:
                axes.errorbar(minutes_grp[num],avg_grp[num],stderr_grp[num],label=filenm[num]+',n='+str(count[num][0]))
        axes.legend(fontsize=10, loc='best')
    fig.canvas.draw()
    pyplot.show()
    return 


def construct_filename(sepvarlist,paramgrp):
    filnm=''
    for sepnum in range(len(sepvarlist)):
        sepvar=sepvarlist[sepnum][0]
        attr=getattr(paramgrp[0],sepvar)
        filnm=filnm+sepvar+str(attr)
	#filnm=filnm+str(attr)
    return filnm

def opto_filename(sepvarlist,paramgrp,llval):
    filnm=''
    for sepnum in range(len(sepvarlist)):
        sepvar=sepvarlist[sepnum][0]
        if sepvar=='genotype':
            sep_phrase='geno'
            attr=getattr(paramgrp[0],sepvar)[0:-1]
        elif sepvar=='cre':
            sep_phrase=sepvar
            attr=getattr(paramgrp[0],'genotype')[-1]
        else:
            attr=getattr(paramgrp[0],sepvar)
            sep_phrase=sepvar
            if sepvar=='lightlevel':
 			#print 'filename', sepvar,attr,
               	sep_phrase='LL'
                if attr> llval:
                    attr=str(llval+1)+'-100'
                else:
                    attr='0-'+str(llval)
            if sepvar=='compound':
                sep_phrase="comp"
                if attr>0:
                    attr='C'
                else:
                    attr='S'
        #filnm=filnm+sep_phrase+str(attr)
        filnm=filnm+str(attr)
	#print filnm
    return filnm
    
def plot_onegroup(dict_group,param_group,group_name):
    fig=pyplot.figure(figsize=(10,15))
    #fig,axes=pyplot.subplots(1,1)
    fig.canvas.set_window_title('popspike vs time for traces of '+group_name)
    ps=[p['popspikenorm'] for p in dict_group]
    pstime=[p['popspikeminutes'] for p in dict_group]
    expername=[p.exper for p in param_group]
    for expname,index in zip(expername,range(len(ps))):
        pyplot.plot(pstime[index],ps[index],label=expname)
    pyplot.legend()
    fig.canvas.draw()
    pyplot.show()

def plot_corr(dict_group,param_group,paramgrpnum,numgroups,filenm,firstpt,lastpt,base_min):
    fig,axes=pyplot.subplots(2,1)
    fig.canvas.set_window_title('normpopspike vs ')
    for gnum in range(numgroups):
        ps=[np.nanmean(p['popspikenorm'][firstpt:lastpt]) for p in dict_group[paramgrpnum[gnum]]]
        epsp=[np.nanmean(p['amp'][0:base_min]) for p in dict_group[paramgrpnum[gnum]]]
        age=[p.age for p in param_group[paramgrpnum[gnum]]]
        for ax,item in enumerate([epsp,age]):
            if np.isnan(ps).any():
                print ("group", filenm[gnum],ps)
            else:
                popt,pcov=optimize.curve_fit(psu.line,item,ps)
                Aopt,Bopt=popt
                Astd,Bstd=np.sqrt(np.diag(pcov))
                labl=str(np.round(Bopt,3))+'+/-'+str(np.round(Bstd,3))
            axes[ax].plot(item,ps,'o',label=filenm[gnum][6:]+labl)
    axes[0].set_xlabel('baseline epsp (mV)')
    axes[1].set_xlabel('age (days)')
    for axis in axes:
        axis.set_ylabel('normPopSpike')
        axis.legend(fontsize=10, loc='best')
    fig.canvas.draw()
    pyplot.show()
