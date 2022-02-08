import numpy as np 
import glob
import pickle
import sys
import PopSpikeAnalClass as psc
from argparse import Namespace
import pop_spike_utilities as psu

ms_per_sec=1000

def find_files(subdir):
	pattern = subdir+'202201*.pickle'
	outfnames = sorted(glob.glob(pattern))
	if not len(outfnames):
		print('************ No files found *********')
	return outfnames
	
def analyze_files(fnames,new_outdir,do_exit):
	#loop over all analyzed exepriments
	slopechange=[]
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
			param_dict = vars(params)

			#3. extract analysis parameters, store all needed parameters in params
			param_dict['base_min']=datadict['anal_params']['baseline_min']
			param_dict['noisethres']=datadict['anal_params']['noisethresh']
			param_dict['decay']=datadict['anal_params']['artdecay']*ms_per_sec
			param_dict['FVwidth']=datadict['anal_params']['FVwidth']*ms_per_sec
			param_dict['artthresh']=datadict['anal_params']['artifactthresh']
			param_dict['art_win']=datadict['anal_params']['artifact_wind']
			param_dict['first_peak_end_frac']=datadict['anal_params']['1st_peak_end']
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
			pop_spike.normalizeAmp()
			pop_spike.summaryMeasure()
			print('PREVIOUS SLOPE for', params.exper,'=',old_slope,'+/-',old_slope_std)
			pop_spike.showSummaryPlot(plot=True)
			if np.abs((old_slope-pop_spike.Bopt)/old_slope)>0.1: #percent change from old slope
				param_dict['old slope']=old_slope
				param_dict['new slope']=round(pop_spike.Bopt,6)
				slopechange.append(param_dict)
			pop_spike.params.exper=pop_spike.params.exper+'_newb'
			pop_spike.params.outputdir=new_outdir
			pop_spike.saveData()  #save data in new directory called, e.g. "fixed baseline"
			#create additional file that lists experiments in which original baseline was good, and now is not
			# or original baseline was bad, and now is good
	return slopechange

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
slopechange=analyze_files(pickle_files,new_outdir,do_exit)
if len(slopechange):
	print(slopechange)
	######### Write text file of experiments that have significant slope changes ####
	variables=['exper','old slope','new slope','sex','age','drug','theta','region','decay','artthresh','FVwidth','first_peak_end_frac','base_min','noisethres']
	header='   '.join(variables)
	f=open('slopechange.txt','w')
	f.write(header+'\n')
	for exper in slopechange:
		newline='   '.join([str(exper[k]) for k in variables])
		f.write(newline+'\n')
	f.close()
else:
	print('NO MAJOR SLOPE CHANGE')
