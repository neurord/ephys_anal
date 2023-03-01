import numpy as np
from matplotlib import pyplot
from matplotlib import gridspec
from scipy import optimize
import pop_spike_utilities as psu

######## This function calculates averages across experiments
def exp_avg(datas, datatype=None,somethreshold=999e6):
    if not datatype:
        alldatavalues=datas
        #print("exp_avg: datatype=none, single array passed")
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

def plot_groups(avg_grp,stderr_grp,minutes_grp,count,filenm,sepvarlist,sepvardict, plot_cols=None):
    #option: 1 or multiple columns, even if two or more seperation variables
    #maxdur=minutes_grp[0][-1]
    keylist=sorted(list(avg_grp.keys()))
    if len(sepvarlist)==1:
        entry1 = list(np.unique([k for k in keylist]))
        entry2=[]
    else:
        entry1 = list(np.unique([k[0] for k in keylist]))
        entry2 = list(np.unique([k[1] for k in keylist]))
    if len(sepvarlist)==1 or plot_cols==1:# or some parameter specified
        numcols=1
    else:
        numcols=len(entry2)
    numrows=len(entry1)
    print('row keys', entry1,'col keys', entry2,'#rows, #cols', numrows,numcols)
    #set up axes for the plots
    pyplot.ion()
    fig=pyplot.figure(figsize=(10,15))
    fig.canvas.set_window_title('Group Average of Normalized PopSpike')
    gs=gridspec.GridSpec(numrows,numcols)
    axes=[]

    for row in range(numrows):
        for col in range(numcols):
            axes.append(fig.add_subplot(gs[row,col]))
            axes[-1].axis([-15,60,.6,1.6])
            axes[-1].axhline(1)
            axes[-1].set_ylabel('normPopSpike')
            axes[-1].set_xlabel('Time (min)')
        
    for grp in avg_grp.keys():  
        if len(sepvarlist)==1:
            axnum=entry1.index(grp)
        elif numcols==1:
            row=entry1.index(grp[0])
            col=0
            axnum=row
        else:
            row=entry1.index(grp[0])
            col=entry2.index(grp[1])
            axnum=row*numcols+col
        axes[axnum].errorbar(minutes_grp[grp],avg_grp[grp],stderr_grp[grp],label=filenm[grp]+',n='+str(np.max(count[grp])))
        axes[axnum].legend(fontsize=10, loc='best')
    return fig

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

def plot_corr(grp_data,samp=1): 
    fig,axes=pyplot.subplots(2,1)
    fig.canvas.set_window_title('normpopspike vs ')
    for grp in grp_data.groups.keys():
        #ps=np.nanmean(grp_data.get_group(grp).popspikenorm[firstpt:lastpt],axis=0)
        ps=[grp_data.get_group(grp).PS_mean.iloc[i][samp] for i in range(len(grp_data.get_group(grp)))]
        epsp=[grp_data.get_group(grp).PS_mean.iloc[i][0] for i in range(len(grp_data.get_group(grp)))]
        age=grp_data.get_group(grp)['age'].values #replace age with list of params?
        for ax,item in enumerate([epsp,age]):
            if np.isnan(ps).any() or len(ps)<3:
                print ("group", grp,ps)
                labl=''
            else:
                popt,pcov=optimize.curve_fit(psu.line,item,ps)
                Aopt,Bopt=popt
                Astd,Bstd=np.sqrt(np.diag(pcov))
                labl=str(np.round(Bopt,3))+'+/-'+str(np.round(Bstd,3))
            axes[ax].plot(item,ps,'o',label=' '.join([str(g) for g in grp])+' ,'+labl)
    axes[0].set_xlabel('baseline epsp (mV)')
    axes[1].set_xlabel('age (days)')
    for axis in axes:
        axis.set_ylabel('normPopSpike')
        axis.legend(fontsize=10, loc='best')
    fig.canvas.draw()
    pyplot.show()
