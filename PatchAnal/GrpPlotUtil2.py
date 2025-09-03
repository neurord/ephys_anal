import numpy as np
np.set_printoptions(legacy='1.25')
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
                print ("plot_corr: Nans in yalues or <3 samples for group:", group,'X:', xvals,'Y:', yvalues)
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

def plot_IVIF_mean(grp,ivif_dict): 
    plot_vars=[a for a in grp.IVIF_variables if a != 'Im']
    fig,axes=pyplot.subplots(len(plot_vars),1,sharex=True)
    axes=fig.axes
    fig.suptitle(', '.join(plot_vars)+' vs '+'current')
    fig.canvas.manager.set_window_title('IV_IF')
    for group in ivif_dict.keys():
        xvals= ivif_dict[group]['Im']#Current
        #print ("group:", group,'Im:',xvals)
        for ax,yvar in enumerate(plot_vars):
            yvalues=ivif_dict[group][yvar]
            ste=ivif_dict[group][yvar+'_ste']
            #print ('yvar',yvar,', Y:',yvalues)
            labl=group_to_word(group)
            axes[ax].errorbar(xvals,yvalues,ste,label=labl,marker='.',capsize=4) #linestyle='',
            #axes[ax].plot(xvals,yvalues,label=labl)
            axes[ax].set_ylabel(yvar)
            axes[ax].set_xlabel('Im')
    for axis in axes:
        axis.legend(fontsize=10, loc='best')
    fig.canvas.draw()
    pyplot.show()

def plot_IVIF(grp,yvars): 
    figs=[]
    plot_variables=[a for a in grp.IVIF_variables if a != 'Im']
    for fnum,yvar in enumerate(plot_variables):
        fig,axes=pyplot.subplots(len(grp.grp_data.groups.keys()),1)
        figs.append(fig)
        axes=fig.axes
        fig.suptitle(yvar+' vs '+'current')
        fig.canvas.manager.set_window_title('IV_IF separate'+str(fnum))
        for ax,group in enumerate(grp.grp_data.groups.keys()):
            labl=group_to_word(group)
            for ii in grp.grp_data.get_group(group).index:
                yvalues=grp.grp_data.get_group(group)[yvar][ii]
                xvals=grp.grp_data.get_group(group).Im[ii]*1e12 #Current
                axes[ax].plot(xvals,yvalues,label=grp.grp_data.get_group(group)['exper'][ii])
                axes[ax].set_ylabel(yvar+' '+labl)
                axes[ax].set_xlabel('Im (pA)')
        for axis in axes:
                axis.legend(fontsize=10, loc='best')
        fig.canvas.draw()
        pyplot.show()
    return figs

def form_clusters(xvars, eps):
    all_xvals=[]
    for i in xvars.index:
        #create tuples: (Im, index)
        all_xvals.append([(xv,i,j) for j,xv in enumerate(xvars[i])])
    #flatten the list
    all_xvals=[ax for one_set in all_xvals for ax in one_set]
    #sort by xval
    sorted_xvals=sorted(all_xvals)
    clusters=[]
    curr_cluster = [sorted_xvals[0]]
    for curr_point,point in zip(sorted_xvals[0:-1],sorted_xvals[1:]):
        if point[0] <= curr_point[0] + eps:
            curr_cluster.append(point)
        else:
            clusters.append(curr_cluster)
            curr_cluster = [point]
    clusters.append(curr_cluster)
    return clusters

def cluster_IVIF (grp,eps=5):     #Use eps=5 or 10 pA
    import scipy.stats
    #use the index (second value of tuple) to guide the clustering.  Might need to add 2nd index of which inject it is
    ivif_dict={a:{} for a in grp.grp_data.groups.keys()}
    yvars=[a for a in grp.IVIF_variables if a != 'Im'] 
    for group in grp.grp_data.groups.keys():
        xvars=grp.grp_data.get_group(group)['Im']*1e12  #convert to pA
        clusters=form_clusters(xvars,eps)
        mean_x=[]
        for clust in clusters:
            mean_x.append(np.mean([c[0] for c in clust]))
        ivif_dict[group]['Im']=mean_x
        for yvar in yvars:
            mean_y=[];ste_y=[]
            yvalues=grp.grp_data.get_group(group)[yvar]
            for clust in clusters:
                mean_y.append(np.nanmean([yvalues[c[1]][c[2]] for c in clust]))
                ste_y.append(scipy.stats.sem([yvalues[c[1]][c[2]] for c in clust],nan_policy='omit'))
            ivif_dict[group][yvar]=mean_y
            ivif_dict[group][yvar+'_ste']=ste_y
    return ivif_dict
            


            


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
                    vars()[iv][each_dict[IDfield]]=each_dict[iv] #add to dictionary
        for iv in indep_var:
            grp.whole_df[iv] = grp.whole_df['ID'].map(vars()[iv])
    print(grp.whole_df[print_vars])
