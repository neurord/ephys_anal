import numpy as np
from matplotlib import pyplot
from matplotlib import gridspec
from scipy import optimize
from PatchAnalPlots import line

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

def plot_groups(avg_grp,stderr_grp,minutes_grp,count,common_filenm,sepvarlist,plot_cols=None):
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
    fig.canvas.manager.set_window_title('Group Average of Normalized PSP for:'+common_filenm)
    gs=gridspec.GridSpec(numrows,numcols)
    axes=[]

    for row in range(numrows):
        for col in range(numcols):
            axes.append(fig.add_subplot(gs[row,col]))
            #axes[-1].axis([-10,40,0,2]) #axis limits! FIXME: Definitely change these
            axes[-1].axhline(1)
            axes[-1].set_ylabel('normPSP')
            axes[-1].set_xlabel('Time (min)') #FIXME: use stim current when plotting IO curve
        
    for grp in avg_grp.keys():  
        if len(sepvarlist)==1:
            axnum=entry1.index(grp)
            lbl=grp
        elif numcols==1:
            row=entry1.index(grp[0])
            col=0
            axnum=row
            lbl='_'.join(grp)+',n='+str(np.max(count[grp]))
        else:
            row=entry1.index(grp[0])
            col=entry2.index(grp[1])
            axnum=row*numcols+col
            lbl='_'.join(grp)+',n='+str(np.max(count[grp]))
        axes[axnum].errorbar(minutes_grp[grp],avg_grp[grp],stderr_grp[grp],label=lbl) #alternative: filenm[grp]
        axes[axnum].legend(fontsize=10, loc='best')
    return fig
    
def plot_onegroup(dict_group,param_group,group_name):
    fig=pyplot.figure(figsize=(10,15))
    #fig,axes=pyplot.subplots(1,1)
    fig.canvas.manager.set_window_title('popspike vs time for traces of '+group_name)
    ps=[p['popspikenorm'] for p in dict_group]
    pstime=[p['popspikeminutes'] for p in dict_group]
    expername=[p.exper for p in param_group]
    for expname,index in zip(expername,range(len(ps))):
        pyplot.plot(pstime[index],ps[index],label=expname)
    pyplot.legend()
    fig.canvas.draw()
    pyplot.show()

def plot_corr(grp,xvar,yvar,samp=1): 
    fig,axes=pyplot.subplots(len(xvar),1)
    axes=fig.axes
    fig.suptitle(yvar+' vs '+' ,'.join(xvar))
    fig.canvas.manager.set_window_title(yvar)
    for ax,item in enumerate(xvar):
        for group in grp.grp_data.groups.keys():
            yvalues=[grp.grp_data.get_group(group)[yvar].values[i][samp] for i in range(len(grp.grp_data.get_group(group)))]
            xvals=grp.grp_data.get_group(group)[item].values #replace age with list of params?
            if np.isnan(yvalues).any() or len(yvalues)<3:
                print ("group:", group,'Time:', xvals,'Y:', yvalues)
                labl=group_to_word(group)
            else:
                popt,pcov=optimize.curve_fit(line,xvals,yvalues)
                Aopt,Bopt=popt
                Astd,Bstd=np.sqrt(np.diag(pcov))
                labl=group_to_word(group)+' ,'+str(np.round(Bopt,3))+'+/-'+str(np.round(Bstd,3))
            axes[ax].plot(xvals,yvalues,'o',label=labl)
            axes[ax].set_xlabel(item)
            axes[ax].set_ylabel(yvar)
    for axis in axes:
        axis.legend(fontsize=10, loc='best')
    fig.canvas.draw()
    pyplot.show()

def plot_IVIF(grp,yvars): 
    import scipy.stats
    fig,axes=pyplot.subplots(len(yvars),1)
    axes=fig.axes
    fig.suptitle(', '.join(yvars)+' vs '+'current')
    fig.canvas.manager.set_window_title('IV_IF')
    for ax,yvar in enumerate(yvars):
        for group in grp.grp_data.groups.keys():
            yvalues=np.nanmean(np.vstack(grp.grp_data.get_group(group)[yvar].values),axis=0)
            ste=scipy.stats.sem(grp.grp_data.get_group(group)[yvar].values,axis=0,nan_policy='omit')
            xvals=grp.grp_data.get_group(group).Im.mean() #Current
            print ("group:", group,'Im:',xvals,'\nY:',yvalues)
            labl=group_to_word(group)+' , n='+str(len(grp.grp_data.get_group(group)))
            axes[ax].errorbar(xvals,yvalues,ste,label=labl)
            axes[ax].set_ylabel(yvar)
            axes[ax].set_xlabel('Im')
    for axis in axes:
        axis.legend(fontsize=10, loc='best')
    fig.canvas.draw()
    pyplot.show()

def group_to_word(grp):
    if isinstance(grp,tuple) or isinstance(grp,list):
        grp_nm='_'.join(grp)
    elif isinstance(grp,str):
        grp_nm=grp
    else:
        print('unknown group structure.  Not tuple or string or list')
        grp_nm=''
    return grp_nm

def read_IDfile(grp,IDfield,indep_var):
    import csv
    print_vars=['exper','ID']
    #initialize one dictonary for each indep var
    for id in indep_var:
        vars()[id]=dict()
        print_vars.append(id)
    with open(grp.subdir+grp.IDfile+'.csv', newline='', encoding="utf-8") as mycsvfile: #Seems that utf-8 encoding doesn't always work
        my_reader = csv.DictReader(mycsvfile)
        column_names = my_reader.fieldnames ##Grabs the column names from the csv file. 
        for each_dict in my_reader:  #each line is a dictionary
            if each_dict[IDfield]: #if ID field is not blank
                for iv in indep_var: #vars()[id] is the dictinary name
                    vars()[iv]['JK-'+each_dict[IDfield]]=each_dict[iv] #add to dictionary
        for iv in indep_var:
            grp.whole_df[iv] = grp.whole_df['ID'].map(vars()[iv])
    print(grp.whole_df[print_vars])
