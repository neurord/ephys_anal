import numpy as np
from matplotlib import pyplot
pyplot.ion()

MV_PER_V=1000
SEC_PER_MIN=60

def line(x,A,B):
     return B*x+A

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
                dt=exp.get_dt(exp.IV_IF_dict[r]) #FIXME - delete this, fix next line
                time = np.arange(0,len(exp.data['Data'][r].__array__()))*dt# exp.IV_dt[r] #time array needed for trace plot
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

