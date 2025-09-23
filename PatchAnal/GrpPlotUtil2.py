import numpy as np
np.set_printoptions(legacy='1.25')
from matplotlib import pyplot
from matplotlib import gridspec
from scipy import optimize
from PatchAnalPlots import line

SEC_PER_MIN =60

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
    from matplotlib.patches import Rectangle
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
    fig=pyplot.figure(figsize=(6*numcols,10))
    fig.canvas.manager.set_window_title('Group Average of Normalized PSP for:'+common_filenm)
    gs=gridspec.GridSpec(numrows,numcols)
    axes=[]

    for row in range(numrows):
        for col in range(numcols):
            axes.append(fig.add_subplot(gs[row,col]))
            #axes[-1].axis([-10,40,0,2]) #axis limits! FIXME: Definitely change these
            axes[-1].axhline(1)
            axes[-1].set_ylabel('normPSP')
            axes[-1].set_xlabel('Time (min)') 
        
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
        axes[axnum].legend(fontsize=10, loc='upper left')
    for axnum in range(len(axes)):
        axes[axnum].add_patch(Rectangle((0,0.5),3,3,alpha=0.2,color='gray'))
    return fig
    
def plot_onegroup(grp,yvars,factors):
    numcols=len(grp.grp_data.groups.keys())
    fig,axes=pyplot.subplots(len(yvars),numcols,figsize=(3*numcols,3.5*len(yvars)))
    axes=fig.axes
    for row,yvar in enumerate(yvars):
        for col,group in enumerate(grp.grp_data.groups.keys()):
            ax=row*numcols+col
            axes[col].set_title(group_to_word(group))
            for i in grp.grp_data.get_group(group).index:
                xvals=grp.grp_data.get_group(group)['psptime'][i]/SEC_PER_MIN  #convert to minutes
                yvals=grp.grp_data.get_group(group)[yvar][i]*factors[row] #convert to reasonable units
                expname=grp.grp_data.get_group(group).exper[i]
                axes[ax].plot(xvals,yvals,label=expname)
            ylim=axes[ax].get_ylim()
            if np.min(ylim)<0:
                axes[ax].set_ylim([round(ylim[0]),0])
            elif np.min(ylim)>0:
                axes[ax].set_ylim([1,round(ylim[1])])
    for col in range(numcols): #only add xlabel and legend to bottom row
        axes[row*numcols+col].set_xlabel('time')
        axes[row*numcols+col].legend()
    for row,yvar in enumerate(yvars):
        axes[row*numcols].set_ylabel(yvar)    
    return fig

def plot_bad(grp):
    reasons={'nan':'too many nans','slope': 'bad slope','num_traces': 'too few follow-up traces','meanpre': 'baseline psp too small'}
    groups=list(np.unique(grp.bad_df.reason))
    fig,axes=pyplot.subplots(len(groups),1,figsize=(5,3.0*len(groups)))
    fig.canvas.manager.set_window_title('Problems')
    axes=fig.axes
    for jj in grp.bad_df.index:
        row=grp.bad_df.iloc[jj]
        ax=groups.index(row.reason)
        axes[ax].set_title(reasons[row.reason]) 
        axes[ax].set_ylabel('normPSP,'+row.reason)
        num_traces=row['num_traces']
        time=row['psptime']/SEC_PER_MIN #convert to minutes
        yvals=row['normPSP'] #convert to reasonable units
        expname=row.exper
        p=axes[ax].plot(time,yvals,label=expname)
        color = p[-1].get_color()
        if row.reason=='nan':
            nan_index=np.argwhere(np.isnan(row['normPSP']))
            axes[ax].plot(time[nan_index],np.ones((len(nan_index))),'o', color=color)
        elif row.reason=='slope':
            #FIXME: this does not plot correct 5 min line when NOT aligning induction with time=0
            axes[ax].plot(time,row['Intercept']/row['meanpre']+row['slope']*(SEC_PER_MIN/row['meanpre'])*time,linestyle=':', color=color)
            axes[ax].plot(time,row['Intercept10']/row['meanpre']+row['slope10']*(SEC_PER_MIN/row['meanpre'])*time,linestyle='-.', color=color)
            axes[ax].set_title('dotted: 5 min slope, dot-dash: 10 min slope')  
    for ax in axes:
        ax.legend(fontsize=8, loc='best')
    axes[-1].set_xlabel('Time (min)')


def plot_corr(grp,xvar,yvar,samp=0): 
    fig,axes=pyplot.subplots(len(xvar),1,figsize=(6,3*len(xvar)))
    axes=fig.axes
    fig.suptitle(yvar+' vs '+' ,'.join(xvar))
    fig.canvas.manager.set_window_title(yvar)
    for ax,x in enumerate(xvar):
        for group in grp.grp_data.groups.keys():
            yvalues=[grp.grp_data.get_group(group)[yvar].values[i][samp] for i in range(len(grp.grp_data.get_group(group)))]
            xvals=grp.grp_data.get_group(group)[x].values 
            plot_line=False
            if np.isnan(yvalues).any() or len(yvalues)<3:
                print ("plot_corr: Nans in yalues or <3 samples for group:", group,'X:', xvals,'Y:', yvalues)
                labl=group_to_word(group)
            else:
                popt,pcov=optimize.curve_fit(line,xvals,yvalues)
                Aopt,Bopt=popt
                Astd,Bstd=np.sqrt(np.diag(pcov))
                labl=group_to_word(group)+' ,'+str(np.round(Bopt,3))+'+/-'+str(np.round(Bstd,3))
                if  (np.abs(Bopt)-2*Bstd>0): #grp.slope_std_factor
                    plot_line=True
            p=axes[ax].plot(xvals,yvalues,'o',label=labl)
            color=p[0].get_color()
            if plot_line:
                axes[ax].plot(xvals,Aopt+Bopt*xvals,linestyle=':',color=color)
            axes[ax].set_xlabel(x)
            axes[ax].set_ylabel(yvar)
    for axis in axes:
        axis.legend(fontsize=10, loc='best')
    fig.canvas.draw()
    pyplot.show()

def plot_IVIF_mean(grp,ivif_dict, plot_vars, x): 
    fig,axes=pyplot.subplots(len(plot_vars),1,sharex=True,figsize=(6,3*len(plot_vars)))
    axes=fig.axes
    fig.suptitle(', '.join(plot_vars)+' vs '+ x)
    fig.canvas.manager.set_window_title(x+' - group mean') #FIXME: customize for IO
    for group in ivif_dict.keys():
        xvals= ivif_dict[group][x] #Current
        #print ("group:", group,'Im:',xvals)
        for ax,yvar in enumerate(plot_vars):
            yvalues=ivif_dict[group][yvar]
            ste=ivif_dict[group][yvar+'_ste']
            #print ('yvar',yvar,', Y:',yvalues)
            labl=group_to_word(group)
            axes[ax].errorbar(xvals,yvalues,ste,label=labl,marker='.',capsize=4) #linestyle='',
            #axes[ax].plot(xvals,yvalues,label=labl)
            axes[ax].set_ylabel(yvar)
            axes[ax].set_xlabel(x)
    for axis in axes:
        axis.legend(fontsize=10, loc='best')
    fig.canvas.draw()
    pyplot.show()

def sort_y_by_x(Y,X):
    newy=[y for _, y in sorted(zip(X, Y))]
    newx=[x for x, _ in sorted(zip(X, Y))]
    return newy, newx

def plot_IVIF(grp,yvars,x,units): 
    numcols=len(grp.grp_data.groups.keys())
    fig,axes=pyplot.subplots(len(yvars),numcols,figsize=(3*numcols,3.5*len(yvars)))
    fig.canvas.manager.set_window_title(x+' - separate') #FIXME: customize for IO
    axes=fig.axes
    for fnum,yvar in enumerate(yvars):
        for col,group in enumerate(grp.grp_data.groups.keys()):
            ax=fnum*numcols+col
            axes[col].set_title(group_to_word(group))
            xlabel=x+' ('+units+')'
            if units=='pA':
                conversion=1e12
            else:
                conversion=1
            for ii in grp.grp_data.get_group(group).index:
                yvalues=grp.grp_data.get_group(group)[yvar][ii]
                xvals=grp.grp_data.get_group(group)[x][ii]*conversion
                yvalues,xvals=sort_y_by_x(yvalues,xvals) #sort x and y values by ordering x values
                axes[ax].plot(xvals,yvalues,label=grp.grp_data.get_group(group)['exper'][ii])
    for col in range(numcols): #only add xlabel and legend to bottom row
        axes[fnum*numcols+col].set_xlabel(xlabel)
        axes[fnum*numcols+col].legend(fontsize=10, loc='best')
    for fnum,yvar in enumerate(yvars):
        axes[fnum*numcols].set_ylabel(yvar)
    pyplot.show()
    return fig

def form_clusters(xvars, eps,conversion):
    all_xvals=[]
    for i in xvars.index:
        #create tuples: (Im, index)
        all_xvals.append([(xv*conversion,i,j) for j,xv in enumerate(xvars[i])])
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

def cluster_IVIF (grp,yvars,x,eps=5,conversion=1):     #Use eps=5 or 10 pA
    import scipy.stats
    #use the index (second value of tuple) to guide the clustering.  Might need to add 2nd index of which inject it is
    ivif_dict={a:{} for a in grp.grp_data.groups.keys()}
    for group in grp.grp_data.groups.keys():
        xvars=grp.grp_data.get_group(group)[x]
        clusters=form_clusters(xvars,eps,conversion)
        mean_x=[]
        for clust in clusters:
            mean_x.append(np.mean([c[0] for c in clust]))
        ivif_dict[group][x]=mean_x
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
