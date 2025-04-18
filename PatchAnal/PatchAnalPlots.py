import numpy as np
from matplotlib import pyplot
pyplot.ion()

MV_PER_V=1000
SEC_PER_MIN=60

def line(x,A,B):
     return B*x+A

def rows_columns(n):
    factors=list([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)
    index=np.argmin([abs(f[0]-f[1]) for f in factors])
    cols,rows=factors[index]
    return cols,rows

def trace_plot(exp): #plot traces, for visual inspection / verify analysis
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
                axes.plot(time[exp.basestartpt],exp.RMP[headstage][trace_num]*MV_PER_V+trace_num*offset, 'kx')
                trace_num+=1
            axes.set_title('o = peak, x =baseline')
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

def summary_plot(exp):
    fig,axes=pyplot.subplots()
    fig.canvas.manager.set_window_title('Summary '+exp.params['exper'])
    for h,color in zip(exp.pspamp.keys(),['b','r']):
        psp_minutes=(exp.psptime[h])/SEC_PER_MIN #start from 0, convert from sec to min
        axes.plot(psp_minutes,exp.normpsp[h],color=color,marker='.',linestyle='None',label=h)
        axes.plot(psp_minutes[0:exp.num_pre],line(exp.psptime[h][0:exp.num_pre],exp.Aopt[h]/exp.meanpre[h],exp.Bopt[h]/exp.meanpre[h]),color=color,alpha=0.5)
        axes.set_xlabel('Time (minutes)')
        axes.set_ylabel ('Normalized PSP amplitude ')
    axes.legend()
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
            artifact=find_peaks(trace,height=50)[0]
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
        axes.set_xlabel(xlbl)
        axes.set_ylabel ('PSP amplitude ')
    axes.legend()
    fig.show()
