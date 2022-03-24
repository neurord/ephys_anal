import numpy as np 
import glob
import pickle
import sys
import PopSpikeAnalClass as psc
import os
import time
import datetime
from argparse import Namespace
import pop_spike_utilities as psu

ms_per_sec=1000
changedate=datetime.date(2021, 5, 12) #date file changed FVwidth and decay time units

def find_files(subdir):
	pattern = subdir+'20201224greendm105g1.pickle'
	outfnames = sorted(glob.glob(pattern))
	if not len(outfnames):
		sys.exit('************ No files found *********')
	return outfnames
	
def analyze_files(fnames,new_outdir,do_exit):
	#loop over all analyzed exepriments
	slopechange=[]
	expererrors=[]
	for outfname in fnames:
		with open(outfname, 'rb') as f:
			#open pickle file
			datadict = pickle.load(f)
			old_slope=round(datadict['trace']['slope'],4)
			old_slope_std=round(datadict['trace']['slope_std'],4)
			#2. extract experimental parameters 
			params = datadict['parameters']

			#re-create commandline and re-parse args
			exper=params.exper
			region=params.region
			sex=params.sex
			drug=params.drug
			age=params.age
			theta=params.theta
			commandline=[exper,sex,str(age),drug,str(theta),region]
			params=psu.parse_args(commandline,do_exit,0)
			print(params.datadir,glob.glob(params.datadir+params.exper+'*'),params.exper)
			param_dict = vars(params)

			#3. extract analysis parameters, store all needed parameters in params
			#if (int(exper[0:8])<20210512) and (exper!='20160801green') and (exper[0:17]!='20160822orangedm1') and (exper!='20161121orange') and (exper!='20170119greendl5') and (exper!='20170120greendm5') and (exper!='20170324greendm105dmso') and (exper!='20170324greendmctrlnitr') and (exper!='20170324orangedmctrldmso') and (exper!='20170329greendm105nitr') and (exper!='20170329greendmctrldmso') and (exper!='20170329orangedm105dmso') and (exper!='20170405greendlctrlsch23390') and (exper!='20180315greendl105') and (exper!='20180612greendm5') and (exper!='20180916orangedm105') and (exper!='20181005orangedm5') and (exper!='20181117greendl5') and (exper!='20190104greendl5') and (exper!='20190110orangedm105') and (exper!='20190223greendm5') and (exper!='20190224orangedl105') and (exper!='20191211orangedm105mpp') and (exper!='20200113greendl5cycl') and (exper!='20200705orangedl5') and (exper[0:8]!='20200707') and (exper!='20200714greendm105') and (exper!='20201031orangedl5') and (exper!='20201220orangedmctrlg1') and (exper!='20201226orangedm105ppt'): #format of parameters changed on this date
			if os.stat(outfname).st_mtime < time.mktime(changedate.timetuple()):
				#if above fails, use os.path.basename(exper[0:8])
				param_dict['decay']=datadict['anal_params']['artdecay']
				param_dict['FVwidth']=datadict['anal_params']['FVwidth']
				print('1st peak end fraction =',datadict['anal_params']['1st_peak_end'])
			elif exper=='20180312orange1':
				param_dict['decay']=datadict['anal_params']['artdecay']/ms_per_sec
				param_dict['FVwidth']=datadict['anal_params']['FVwidth']/ms_per_sec
				param_dict['first_peak_end_frac']=0.625
			else:
				param_dict['decay']=datadict['anal_params']['artdecay']*ms_per_sec
				param_dict['FVwidth']=datadict['anal_params']['FVwidth']*ms_per_sec
			if datadict['anal_params']['1st_peak_end']<0.5: #this number is arbitrary. maybe use 0.625?
				param_dict['first_peak_end_frac']=0.625
			else:
				param_dict['first_peak_end_frac']=datadict['anal_params']['1st_peak_end']
			param_dict['base_min']=datadict['anal_params']['baseline_min']
			param_dict['noisethres']=datadict['anal_params']['noisethresh']
			param_dict['artthresh']=datadict['anal_params']['artifactthresh']
			param_dict['art_win']=datadict['anal_params']['artifact_wind']
			param_dict['induction_gap_time']=datadict['anal_params']['induction_gap']
			param_dict['baselinepoints']=datadict['anal_params']['baselinepts']

			params = Namespace(**param_dict)
			print('re-analyzing experiment:', params.exper)
			#
			#4. call PopSpikeAnalClass with the parameters, but do not generate plots
			pop_spike = psc.PopSpike(params) 
			pop_spike.read_datafile(question=False)
			#optional: depending on which baseline method we use
			#pop_spike.baselinepoints=int(0.004/pop_spike.dt) #4 msec used by labview
			
			pop_spike.find_induction_gap()
			pop_spike.findPopSpike()
			if len(pop_spike.problemtraces['tracenum'])>5:
				param_dict['error_number']=len(pop_spike.problemtraces['tracenum'])
				expererrors.append(param_dict)
			pop_spike.normalizeAmp()
			pop_spike.summaryMeasure()
			print('PREVIOUS SLOPE for', params.exper,'=',old_slope,'+/-',old_slope_std)
			pop_spike.showSummaryPlot(plot=True)
			if (np.abs(old_slope)>0.01 and np.abs(pop_spike.Bopt)<0.01) or (np.abs(old_slope)<0.01 and np.abs(pop_spike.Bopt)>0.01): #percent change from old slope
				param_dict['old slope']=old_slope
				param_dict['new slope']=round(pop_spike.Bopt,6)
				slopechange.append(param_dict)
			pop_spike.params.exper=pop_spike.params.exper+'_newb'
			pop_spike.params.outputdir=new_outdir
			pop_spike.saveData()  #save data in new directory called, e.g. "fixed baseline"
			#create additional file that lists experiments in which original baseline was good, and now is not
			# or original baseline was bad, and now is good
	return expererrors,slopechange

def writefile(variables,datatoprint,filename):
	header='   '.join(variables)
	f=open(filename,'w')
	f.write(header+'\n')
	for exper in datatoprint:
		newline='   '.join([str(exper[k]) for k in variables])
		f.write(newline+'\n')
	f.close()
	return

ARGS='C:\\Users\\vlewitus\\Documents\\Python_Scripts\\Pickle\\ C:\\Users\\vlewitus\\Documents\\Python_Scripts\\Pickle_new_baseline\\'			
if __name__=='__main__':
    try:
        commandline = ARGS.split() # ARGS="pickle_dir new_outputdir "
        do_exit = False
    except NameError: 
        commandline = sys.argv[1:]
        do_exit = True    

pickledir=commandline[0]
new_outdir=commandline[1]
pickle_files=find_files(pickledir)
expererrors,slopechange=analyze_files(pickle_files,new_outdir,do_exit)

######### Write text file of experiments that have significant slope changes ####
if len(slopechange):
	print(slopechange)
	variables=['exper','old slope','new slope','sex','age','drug','theta','region','decay','artthresh','FVwidth','first_peak_end_frac','base_min','noisethres']
	writefile(variables,slopechange,'slopechange.txt')
else:
	print('NO MAJOR SLOPE CHANGE')
######### Write text file of experiments that have many errors ####
if len(expererrors):
	variables=['exper','error_number','sex','age','drug','theta','region','decay','artthresh','FVwidth','first_peak_end_frac','base_min','noisethres']
	writefile(variables,expererrors,'expererrors.txt')