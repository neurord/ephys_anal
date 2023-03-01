import numpy as np
import sys  
import pandas as pd 
import glob
import pickle
import os
print (os.getcwd())
from matplotlib import pyplot
import pop_spike_utilities as psu
import GrpPlotUtil2 as grp_utl
from GrpAvgPopSpikeClass2 import GrpPopSpike as gps

def read_data(grp_avg):
    DATAS = []
    PARAMS = pd.DataFrame()
    ANAL = pd.DataFrame()
    for outfname in grp_avg.outfnames:
        with open(outfname, 'rb') as f:
            datadict = pickle.load(f)
        exper_param = datadict['parameters']
        numnan=sum(np.isnan(datadict['trace']['popspikenorm']))
        ############ Select subset of files based on user specified criteria ##########
        ignore1 = ((grp_avg.params.sex and grp_avg.params.sex != exper_param.sex) or 
                        (grp_avg.params.age is not None and grp_avg.params.age >= exper_param.age) or
                        (grp_avg.params.maxage is not None and grp_avg.params.maxage <= exper_param.age) or
                        (grp_avg.params.drug and grp_avg.params.drug != exper_param.drug) or
                        (grp_avg.params.region and grp_avg.params.region != exper_param.region) or
                        (grp_avg.params.theta and grp_avg.params.theta != exper_param.theta))
        if ignore1:
            if grp_avg.print_info:
               print ("***** ignoring:", datadict['parameters'])
            next 
        else: #data meets all criteria
            if np.isnan(datadict['trace']['popspikenorm']).any(): #inform about nans, even if using data
                print ("########## np.nan detected", numnan,'times in',exper_param.exper, datadict['anal_params'])
            DATAS.append(datadict['trace'])
            PARAMS = PARAMS.append(vars(datadict['parameters']), ignore_index = True)
            ANAL=ANAL.append(datadict['anal_params'],ignore_index=True)
        
        df2 = pd.DataFrame(DATAS, dtype=float)
        grp_avg.whole_df=pd.concat([PARAMS,ANAL, df2], axis = 1)
        print (len(grp_avg.whole_df),'experiments met criteria.')

def write_traces(gps,current_dict):
    for gnum in gps.grp_data.groups.keys():
        print(gps.grp_data.groups[gnum])
        f=open(gps.filenm[gnum]+"_IO.txt",'w')                      
        outputdata=np.column_stack((current_dict[gnum], gps.samples[gnum],
                            gps.avg_PS[gnum], gps.stderr_PS[gnum]))
        header="current "+gps.filenm[gnum]+"count "+gps.filenm[gnum]+"PS_AVG "+gps.filenm[gnum]+"PS_SE "
        f.write(header+"\n")
        np.savetxt(f,outputdata, fmt='%7.5f',delimiter='   ') #'%7.4f' = format is float with 7 characters, 4 after decimal
        print(f,gps.filenm[gnum],header)
        f.close()

#ARGS = "-outputdir C:\\Users\\kblackw1\\Documents\\Python_Scripts\\ephys_anal\\PopSpike2\\ -sepvarlist sex region theta -plot_ctrl 000"      
if __name__ =='__main__':        

	try:
	    commandline = ARGS.split() #in python: define space-separated ARGS string
	    do_exit = False
	except NameError: #if you are not in python, read in filename and other parameters
	    commandline = sys.argv[1:]
	    do_exit = True
	
	params=psu.parse_args(commandline,do_exit,1)
	print(params, commandline)            
	Grp_PS=gps(params)
	Grp_PS.pattern = Grp_PS.subdir+'*IO.pickle'
	Grp_PS.outfnames = sorted(glob.glob(Grp_PS.pattern))
	print('pattern',Grp_PS.pattern,'files',Grp_PS.outfnames)
	if len(Grp_PS.outfnames):
	    #Grp_PS.read_data() #to use this, need PopSpikeIO.py to save slope
	    read_data(Grp_PS)
	    newcolumnname='drug_combined'
	    Grp_PS.group_data(Grp_PS.sepvarlist,newcolumnname) 
	    current=Grp_PS.whole_df.current.mean()
	    current_dict={k:current for k in Grp_PS.avg_PS.keys()}
	    write_traces(Grp_PS,current_dict)
	    #
	    ###### Graph ###############
	    if int(params.plot_ctrl[0])>0:
	        plot_cols=int(params.plot_ctrl[0])
	    else:
	        plot_cols=None
	    max_ps=np.nanmax([np.nanmax(avg_ps) for avg_ps in Grp_PS.avg_PS.values()])
	    fig=grp_utl.plot_groups(Grp_PS.avg_PS,Grp_PS.stderr_PS,current_dict,Grp_PS.samples,Grp_PS.filenm,Grp_PS.sepvarlist,Grp_PS.sepvardict,plot_cols)
	    axes=fig.axes
	    for ax in axes:
	        ax.axis([0,np.max(current),0,np.ceil(max_ps)])
	        ax.set_xlabel('current (mA)')
	        ax.set_ylabel('PopSpike (mV)')
	    