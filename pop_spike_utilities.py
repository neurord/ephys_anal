import numpy as np
import argparse
import pickle

def parse_args(commandline,do_exit,flag):
    #Descriptive variables used in subsequent program to select similar data to average:
    parser = argparse.ArgumentParser()
    exp=['exper','--exper']
    sex=['sex','--sex']
    age=['age','--age']
    drug=['drug','--drug']
    theta=['theta','--theta']
    region=['region','--region']
    parser.add_argument(exp[flag], type=str, help = 'give exp name, for example: 031011_5Or')
    parser.add_argument('--no-graphs', '-g', dest='graphs', default=True, action='store_false') # -g optional and not position-defined (its absence OK too)
    parser.add_argument(sex[flag], type=str, choices=["M","F", "Fe","Fx"],help="male M or female F or Fe")
    parser.add_argument(age[flag], type=int, help="animal age in days")
    parser.add_argument(drug[flag], type=str, help="what drugs were in the ACFS")
    parser.add_argument(theta[flag], type=float, help="theta frequency (e.g. 5, 8, 10.5), enter 0 for no stim ctrl")   
    parser.add_argument(region[flag], type=str, choices=["DM", "DL"],help="dorsomedial: DM or dorsolateral: DL")
    
    # parser.add_argument('-sample_times', '--list', default = [30,60,90,100], action='append', help = "Enter to modify sample_times (default = [30,60,90,100]),", required = True) 
    #to plot PopSpike vs time for each experiment in a group, make 2nd number 1
	#to plot  correlation between LTP (at summarytime) and age or baseline epsp, make 3rd number 1
    #to plot a single column, make 1st number 1
    parser.add_argument('-plot_ctrl',type=str,default='000') 
    parser.add_argument("-samp_time", nargs="+",default=[30,60,90,120],)   
    parser.add_argument("-sepvarlist", nargs="+",default=['sex','region'],help='list of separation variables for grouping data')
    parser.add_argument('-decay', type = float,  default =1.1, help="artifact decay time (default: 1.7 msec)")
    parser.add_argument('-artthresh',  type = float, default = 1.5, help="artifact threshold (default: 1.5 mV")
    parser.add_argument('-FVwidth',  type = float, default = 1, help="Fiber volley width (default: 1 msec")
    parser.add_argument('-samp_win',  type = int, default =5, help="Average over this many minutes for time samples (default: 5)")
    parser.add_argument('-file_end',  type = str, default = ".Lvm", help="filename ending (default: .Lvm)")
    parser.add_argument('-first_peak_end_frac',  type = float, default = 0.625, help="first_peak_end_fraction (default: 0.625)")
    parser.add_argument('-art_win', type = float, default = 0.5, help="artifact window (default: 0.5)")
    parser.add_argument('-base_min', type = int, default = 15, help="number of minutes for stable baseline (default: 15)")
    parser.add_argument('-noisethres', type = float, default = 1.5, help="Noise threshold (default: 1.5)")
    parser.add_argument('-slope_std', type = int, default = 2, help="Enter to modify slope_std_factor (default: 2)")
    parser.add_argument('-datadir',type = str, default = "G:\\FIELDS\\ValerieData\\",   help="Directory containing raw data")
    parser.add_argument('-outputdir',  type = str, default = "C:\\Users\\vlewitus\\Documents\\Python_Scripts\\Pickle_new_baseline\\", help="Directory for pickle files")
    parser.add_argument('-baselinestart',  type = float, help="time of baseline start in minutes since exper start")
    parser.add_argument('-criticaltimes', type = float, nargs="+", help="time of drug start or changing stimulation current to compensate drug effect")
    if flag:
        parser.add_argument('--maxage', type=int, help="maximum age of animals")

    try:
        args = parser.parse_args(commandline) # maps arguments (commandline) to choices, and checks for validity of choices.
    #if arguments are mapped incorrectly, python wants to exit, but the next line says "don't", instead check whether we are in python (do_exit=False) then don't exit, just give us a warning
    except SystemExit:
        if do_exit:
            raise # raise the exception above (SystemExit) b/c none specified here
        else:
            raise ValueError('invalid ARGS')
    
    #print experiment characteristics to double check unconstrained entries
    print ("exper={}".format(args.exper))
    print ("sex={}".format(args.sex))
    print ("age={}".format(args.age))
    print ("drug={}".format(args.drug))
    print ("theta={}".format(args.theta))
    print ("region={}".format(args.region))
    return args

def line(x,A,B):
    return B*x+A
        
def read_labview(filename):
    with open(filename,'r') as f:
        for line in f:
            if "Samples" in line:
                break
    timepoints_per_trace=int(line.split()[1])
    data = np.loadtxt(filename, skiprows=24) 
    #first column of datafile is time, second column is Vm
    wholetime=data[:,0]
    tempVm=data[:,1]
    dt=wholetime[1]-wholetime[0]
    xtime=wholetime[0:timepoints_per_trace]
    #reshape the data to put each trace in separate column
    datalength=np.shape(tempVm)[0]
    numtraces=int(datalength/timepoints_per_trace)
    Vm_traces=np.reshape(tempVm, (numtraces, -1))
    #extract the start time for each trace
    tracetime=np.zeros(numtraces)
    for i,j in enumerate(range(0, datalength, timepoints_per_trace)):
        tracetime[i]=wholetime[j]
    return Vm_traces,tracetime,dt,xtime