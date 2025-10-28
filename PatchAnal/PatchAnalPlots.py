import numpy as np
from matplotlib import pyplot
pyplot.ion()

MV_PER_V=1000
SEC_PER_MIN=60
MS_PER_S=1e3
NA_PER_AMP=1e9
artifact_height=30 #mV - delta V between two sample points

def line(x,A,B):
     return B*x+A

def rows_columns(n):
    factors=list([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)
    index=np.argmin([abs(f[0]-f[1]) for f in factors])
    cols,rows=factors[index]
    return cols,rows

def trace_plot(exp,ss_color='gray'): #plot traces, for visual inspection / verify analysis
    colors=pyplot.get_cmap('plasma')
    partial_scale=0.9 #avoid the very light 
    offset=1 #in mV, to visualize multiple traces
    colinc={h:(len(colors.colors)-1)/(exp.num_traces[h]-1) if exp.num_traces[h]>1 else 1 for h in exp.headstages }
    for headstage in exp.Vm_psp.keys():  #intersection of headstages requested and those available
        fig,axes=pyplot.subplots()
        fig.canvas.manager.set_window_title('traces '+exp.params['exper']+' '+headstage)
        trace_num=0
        for r in exp.Vm_psp[headstage].keys():
            for num in range(np.shape(exp.data['Data'][r].__array__())[-1]):  #num indexes arrays - resets to zero for each routine
                trace=exp.data['Data'][r].__array__()[:,num]*MV_PER_V 
                time = np.arange(0,len(exp.data['Data'][r].__array__()))*exp.dt #time array needed for trace plot
                plot_pt=int(exp.plotstart/exp.dt) #start plotting data from here
                labl=r.split('_')[0]+'_'+str(num)
                color_index=int(trace_num*colinc[headstage]*partial_scale)
                mycolor=colors.colors[color_index]
                axes.plot(time[plot_pt:],trace[plot_pt:]+trace_num*offset,label=labl,color=mycolor)
                ###### label the peak and baseline
                yval=(exp.pspamp[headstage][trace_num]+exp.RMP[headstage][trace_num])*MV_PER_V +trace_num*offset
                xval=time[int(exp.peaktime[headstage][trace_num]/exp.dt)]
                axes.plot(xval,yval, 'ko')
                hyper_Vm=-exp.Raccess[headstage][trace_num]*np.mean(exp.Iaccess)+exp.RMP[headstage][trace_num]
                if ss_color:
                    base_size=len(time[exp.hyper_endpt:exp.base_endpt])
                    hyper_size=len(time[exp.hyperstartpt:exp.hyper_endpt])
                    hyper_Vm=-exp.Raccess[headstage][trace_num]*np.mean(exp.Iaccess)+exp.RMP[headstage][trace_num]
                    axes.plot(time[exp.hyper_endpt:exp.base_endpt],np.full(base_size,exp.RMP[headstage][trace_num]*MV_PER_V+trace_num*offset),ss_color)
                    axes.plot(time[exp.hyperstartpt:exp.hyper_endpt],np.full(hyper_size,hyper_Vm*MV_PER_V+trace_num*offset),ss_color)
                else:
                    axes.plot(time[exp.basestartpt],exp.RMP[headstage][trace_num]*MV_PER_V+trace_num*offset, 'x',color='gray')
                    axes.plot(time[exp.hyperstartpt],hyper_Vm*MV_PER_V+trace_num*offset,'x',color='gray')
                axes.plot(exp.anal_params['dV_2ms_time'][1]+exp.dt,exp.dV_2ms[headstage][trace_num]*MV_PER_V+trace[int(exp.anal_params['dV_2ms_time'][0]/exp.dt)]+trace_num*offset,'r+')
                axes.plot(exp.anal_params['dV_2ms_time'][0]+exp.dt,trace[int(exp.anal_params['dV_2ms_time'][0]/exp.dt)]+trace_num*offset,'r+')
                trace_num+=1
            axes.set_title('o = peak, x =baseline, +=2 ms pulse')
            axes.set_xlabel('Time (s)')
            axes.set_ylabel('Vm (mV)')
            axes.legend()

def IVIF_plot(exp): #plot traces, for visual inspection / verify analysis
    colors=pyplot.get_cmap('plasma')
    partial_scale=0.9 #avoid the very light 
    for headstage in exp.Vm_IV_IF.keys():  #intersection of headstages requested and those available
        fig,axes=pyplot.subplots(nrows=len(exp.Vm_IV_IF[headstage].keys()),ncols=1)
        fig.canvas.manager.set_window_title('traces '+exp.params['exper']+' '+headstage)
        for i,r in enumerate(exp.Vm_IV_IF[headstage].keys()):
            ylbl='_'.join([r.split('_')[0]]+r.split('_')[2:])
            num_traces=np.shape(exp.data['Data'][r].__array__())[-1]
            colinc=(len(colors.colors)-1)/(num_traces-1) #IF and IV have different numbers of traces
            for num in range(num_traces):  #num indexes arrays - resets to zero for each routine
                trace=exp.data['Data'][r].__array__()[:,num]*MV_PER_V 
                time = np.arange(0,len(exp.data['Data'][r].__array__()))*exp.IV_dt[r] #time array needed for trace plot
                mycolor=colors.colors[int(num*colinc*partial_scale)]
                axes[i].plot(time,trace,label=str(num),color=mycolor)
                axes[i].set_ylabel(ylbl)
                #axes[i].legend()
        axes[i].set_xlabel('Time (s)')
        fig.suptitle('headstage '+headstage)

def twin_plot(xval,yval,col,ax,lbl='',unit=''):
    axtwin = ax.twinx()
    axtwin.scatter(xval*NA_PER_AMP,yval,color=col,label=lbl)
    axtwin.spines['right'].set_visible(True)
    axtwin.spines['right'].set_color(col)
    axtwin.spines['left'].set_visible(False) 
    axtwin.tick_params('y',labelcolor=col,left=False,labelleft=False)
    axtwin.set_ylabel(lbl+' ('+unit+')',color=col)
    return axtwin

def IVIF_measures(exp):
    IVcolor={'H1':'b','H2':'r'}
    figIV,axIV=pyplot.subplots()
    figIV.suptitle('IV curve '+exp.params['exper'])
    axIV.set_xlabel('current (nA)')
    axIV.set_ylabel('Vm (mV)')
    for h,IV_IF in exp.IV_IF.items():
        for r,IVIFset in IV_IF.items():
            if np.max(IVIFset.num_spikes)>0:
                figIF,axes=pyplot.subplots(3,1)
                figIF.suptitle('IF measures '+exp.params['exper'])
                for ax,meas1,meas2,unit in zip(axes,['latency','APwidth','AHP_time'],['APfreq','APheight','AHP_amp'],['Hz','mV','mV']):
                    ax.scatter(IVIFset.Im*NA_PER_AMP,IVIFset[meas1]*MS_PER_S,color='k')
                    ax.set_ylabel(meas1+' (ms)')
                    ax0twin=twin_plot(IVIFset.Im,IVIFset[meas2],'red',ax,meas2,unit)
                axes[-1].set_xlabel('current (nA)')
            else:
                axIV.scatter(IVIFset.Im*NA_PER_AMP,IVIFset.Vm*MV_PER_V,label=h,color=IVcolor[h])
    return

def summary_plot(exp):
    fig,axes=pyplot.subplots(2,1)
    axes=fig.axes
    fig.canvas.manager.set_window_title('Summary '+exp.params['exper'])
    for h,color in zip(exp.pspamp.keys(),['b','r']):
        psp_minutes=(exp.psptime[h])/SEC_PER_MIN #start from 0, convert from sec to min
        axes[0].plot(psp_minutes,exp.normpsp[h],color=color,marker='.',linestyle='None',label=h)
        axes[0].plot(psp_minutes[0:exp.num_pre],line(exp.psptime[h][0:exp.num_pre],exp.Aopt[h]/exp.meanpre[h],exp.Bopt[h]/exp.meanpre[h]),color=color,alpha=0.5)
        nan_list(psp_minutes,exp.normpsp[h],1.0,axes[0],h)
        axes[0].set_ylabel ('Normalized PSP amplitude ')
        axes[0].legend()
        #Now plot RMP and Raccess
        axes[1].plot(psp_minutes,-exp.RMP[h]*MV_PER_V,color=color,marker='.',linestyle='None',label=h+' -RMP (mV)')
        axes[1].plot(psp_minutes,exp.Raccess[h]/1e6,color=color,marker='x',linestyle='None',label=h+' Raccess (MOhm)')
        axes[1].plot(psp_minutes,-exp.dV_2ms[h]*MV_PER_V,color=color,marker='+',linestyle='None',label=h+' dV_2ms (mV)')
        axes[-1].set_xlabel('Time (minutes)')
        axes[1].legend()
    fig.show()

def induction_plot(exp,stim_time=[]): #plot traces, for visual inspection / verify analysis
    def induct_summary(ylabel,attrib,AP=None):
        fig,axes=pyplot.subplots()
        for i,r in enumerate(exp.num_spikes.keys()):
            color_index=int(i*colinc*partial_scale)
            mycolor=colors.colors[color_index]
            index=list(exp.__getattribute__(attrib)[r].keys())
            values=list(exp.__getattribute__(attrib)[r].values())
            if len(values):
                if isinstance(values[0],np.recarray) and AP:
                    values=[x.__getattribute__(AP)[0]for x in values] #could specify WHiCH AP instead of [0]
            axes.scatter(index,values,label=r.split('_')[0],s=(20-2*i)**2,color=mycolor)
        axes.legend()
        axes.set_xlabel('burst number')
        axes.set_ylabel(ylabel)
        return fig

    from scipy.signal import find_peaks
    colors=pyplot.get_cmap('plasma')
    partial_scale=0.9 #avoid the very light 
    offset=2 #in mV, to visualize multiple traces
    n=len(exp.spikes.keys())
    cols,rows=rows_columns(n)
    fig,axes=pyplot.subplots(rows,cols)
    axes=fig.axes
    fig.canvas.manager.set_window_title('induction '+exp.params['exper'])
    for i,r in enumerate(exp.spikes.keys()):
        colinc=(len(colors.colors)-1)/(np.shape(exp.data['Data'][r].__array__())[-1]-1)
        for num in range(np.shape(exp.data['Data'][r].__array__())[-1]):  #num indexes arrays 
            trace=exp.data['Data'][r].__array__()[:,num]*MV_PER_V 
            time = np.arange(0,len(exp.data['Data'][r].__array__()))*exp.induct_dt #time array needed for trace plot FIXME: wrong dt
            labl=str(num)
            color_index=int(num*colinc*partial_scale)
            mycolor=colors.colors[color_index]
            axes[i].plot(time,trace+num*offset,label=labl,color=mycolor)
            ###### label the stim artifacts  and stim time
            stim_pts=[int(st/exp.induct_dt) for st in stim_time]
            yval=(exp.data['Data'][r].__array__()[stim_pts,num])*MV_PER_V +num*offset
            xval=time[stim_pts]
            axes[i].plot(xval,yval, 'ko')
            artifact=find_peaks(trace,height=artifact_height)[0]
            yval=trace[artifact]+num*offset
            xval=time[artifact]
            axes[i].plot(xval,yval,'g*')
            #label the action potentials
            if num in exp.spikes[r].keys():
                yval=exp.spikes[r][num].APpeak*MV_PER_V +num*offset
                xval=exp.spikes[r][num].APtime
                axes[i].plot(xval,yval,'P', color='gray')
        axes[i].set_xlabel('Time (s)')
        axes[i].set_ylabel(r.split('_')[0]+' '+'Vm (mV)')
    axes[0].legend()
    axes[0].set_title('o=~stim time, *=artifact, +=spike')
    fig2=induct_summary('num spikes','num_spikes')
    fig3=induct_summary('1st spike time (s)','spikes',AP='APtime')

def nan_list(xvals,yvals,nanval,ax,h):
    nans=np.where(np.isnan(yvals))[0]
    if len(nans):
        newx=[xvals[x] for x in nans]
        newy=[nanval for x in nans]
        ax.plot(newx,newy,color='black',linestyle='None',marker='*',label='nan'+h)
    return 

def IO_plot(exp):
    fig,axes=pyplot.subplots()
    fig.canvas.manager.set_window_title('IO curve '+exp.params['exper'])
    for h,color in zip(exp.IOamp.keys(),['b','r']):
        if len(exp.IOamp[h])==len(exp.IOrange):
            xvals=exp.IOrange
            xlbl='current'
        else:
            xvals=np.arange(len(exp.IOamp[h]))
            xlbl='index'
        axes.plot(xvals,exp.IOamp[h],color=color,marker='.',linestyle='None',label=h)
        nan_list(xvals,exp.IOamp[h],-.001,axes,h)
        axes.set_xlabel(xlbl)
        axes.set_ylabel ('PSP amplitude ')
    axes.legend()
    fig.show()

def plot_hist(measures,dfset,exp):
    from matplotlib import pyplot as plt
    plt.ion()
    n=len(measures)
    cols,rows=rows_columns(n)
    fig,axes=plt.subplots(rows,cols)
    axes=fig.axes
    for key in dfset.keys():
        for ax,m in enumerate(measures):
            minv=np.percentile(dfset[key][m].dropna().values,1)
            maxv=np.percentile(dfset[key][m].dropna().values,99)
            axes[ax].hist(dfset[key][m].values,range=(minv,maxv),bins=exp.bins)
            if m in exp.units:
                label=m+'*'+exp.units[m]
            else:
                label=m
            axes[ax].set_xlabel(label)
    fig.canvas.manager.set_window_title('Histogram '+exp.params['exper'])
    fig.tight_layout()
    return fig

def plot_peaks(fig,axes, time,amp,color,symbol=None):
    for i,(tm,am) in enumerate(zip(time,amp)):
        axes[i].scatter(tm,am,marker=symbol,c=color) #use text marker in plot_traces

def plot_pairs(dfset,pairs):
    from matplotlib import pyplot as plt
    plt.ion()
    cols,rows=rows_columns(len(pairs))
    fig,axes=plt.subplots(rows,cols)
    axes=fig.axes
    for ax,p in enumerate(pairs):
        for key in dfset.keys():
            corr=dfset[key][p[0]].corr(dfset[key][p[1]])
            axes[ax].scatter(dfset[key][p[0]],dfset[key][p[1]],label=key+' '+str(round(corr,4)))
            axes[ax].set_xlabel(p[0])
            axes[ax].set_ylabel(p[1])
            axes[ax].legend()            
        
