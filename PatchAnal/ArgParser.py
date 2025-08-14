import argparse

def ArgParserPatch(commandline,do_exit,flag=0):
    parser = argparse.ArgumentParser()
    exp=['experiment','--exper']
    parser.add_argument(exp[flag], type=str, help = 'give name of h5 file, for example: ')
    parser.add_argument('-headstages', type=str, nargs="+", help='which headstage to analyze (default=both)', default=['H1','H2']) #FIXME - make this required for PatchAna?
    parser.add_argument("-celltype", type=str, nargs="+", choices=["D1-SPN", "FSI", "D2-SPN",'ChI','other'],help="D1-SPN, D2-SPN, FSI or other cell type, one for each headstage") #FIXME - make this required for PatchAnal?
    parser.add_argument('-digstim',type=float,help='value of dig stim for theta') #Value of digi stimulation. FIXME - make this required
    parser.add_argument('-datadir',type = str, default = "C:\\Users\\klblackwell\\Documents\\Python\\ephys_anal\\PatchAnal\\Theta\\",   help="Directory containing raw data")
    parser.add_argument('-outputdir',  type = str, default = "C:\\Users\\klblackwell\\Documents\\Python\\ephys_anal\\PatchAnal\\Theta\\", help="Directory for npz files")
    parser.add_argument('-filetype',  type = str, default = ".h5", help="filetype of experiment file")
    parser.add_argument('--no-graphs', '-g', dest='graphs', default=True, action='store_false') # -g optional and not position-defined (its absence OK too)
    parser.add_argument('-plotstart',type=float,default=0,help='time in sec when plots will start')
    ### next set = meta data - read from meta data in file?
    #parser.add_argument("-bathtemp", type=str, choices=["heat","RT"],help="heat or RT", default='RT')
    #parser.add_argument("-Rtip", type=float,help="pipette tip resistance") #Read from logfile
    parser.add_argument("-ID", type=str, help="animal ID") #READ FROM NOTEBOOK file 
    parser.add_argument("-sex", type=str, help="sex/estrus/pellet, e.g. M, F, FE, OVX, OVX+E2") #READ FROM NOTEBOOK file for PatchAnal, no default for GrpAvgPatchClass 
    parser.add_argument("-age", type=int, help="animal age in days") #READ FROM NOTEBOOK file  for PatchAnal, no default for GrpAvgPatchClass 
    parser.add_argument("-drug", type=str,help="bath solution") #READ FROM NOTEBOOK file for PatchAnal, no default for GrpAvgPatchClass 
    parser.add_argument("-region", type=str, choices=["DMS", "DLS"], help="which striatal region") # READ FROM NOTEBOOK file for PatchAnal, no default for GrpAvgPatchClass 
    parser.add_argument("-genotype", type=str, choices=["tdTomato", "wt"],help="tdTomato or wt", default='tdTomato') #eventually read from metadata
    parser.add_argument('-slope_std', type = int, default = 2, help="baseline slope_std_factor criteria (default: 2)")
    parser.add_argument('-window', type=int, help='number of points to average for peak Vm', default=5) #also used for time samples
    #next set refers to the IV, IF curves and the chirp - likely these are not needed as can be read from the hdf5 file
    #parser.add_argument('-chirp', nargs= '+', default=[4,5], help='series number for chirps')
    #parser.add_argument('-IF', nargs= '+', default = [1, 40e-12, 20e-12], help='IF series number, starting point, increment')
    #parser.add_argument('-IV', nargs='+', default = [2, -450e-12, 50e-12], help='IV series number, starting point, increment')
    #parser.add_argument('-injstart', type=float, help="injection start time", default= 0.2)
    #parser.add_argument('-injstop', type = float, help="injection stopping time", default = 0.6)
    if not flag:
        parser.add_argument('-IOrange',type=float,nargs='+',help='values of dig stim for IOtest, space separated',default=[0.02, 0.12, 0.22, 0.32, 0.42, 0.52, 0.62, 1.2, 2.2, 3.2])
        parser.add_argument('-PSPstart', type=float, help = 'earliest time a PSP could be detected', default=0.927)
        parser.add_argument('-basestart', type=float, help='time PRIOR to event for assessing baseline membrane potential', default=0.12)
        parser.add_argument('-base_dur', type=float, help='duration of baseline period, make smaller than basestart', default=0.10) # <= basestart
        parser.add_argument('-ss_dur', type=float, help='calculate steady state Vm using last X sec of pulse', default=0.050)
        parser.add_argument('-APthresh', type=float,help='minimum amplitude to be considered spike',default=0)
        parser.add_argument('-thresheight', type=float,help='minimum spike height to be considered spike',default=0.05)
        parser.add_argument('-refract', type=float,help='minimum time between AP',default=0.004)
        parser.add_argument('-max_risetime',type=float,help='AP rejected if risetime slower than this',default=0.004)
        parser.add_argument('-min_risetime',type=float,help='AP rejected if risetime faster than this',default=0.0003)
        parser.add_argument('-threshval',type=float,help='Vm at earliest point that rise exceeds 2% of max risetime', default=0.02) #with 5%, bigger AHP, smalle spike
        parser.add_argument('-PSP_interval',type=int,help='time between PSP baseline traces, if no notebook file', default=30) #with 5%, bigger AHP, smalle spike
        parser.add_argument("-decay", type = float, help='stimulation artifact decay time, in sec', default = 0.002) #includes time for AP if one occurs
        parser.add_argument("-induction", type=str, choices=["ThetaBurst", "20Hz", 'None'], help='name of induction protocol', default = 'ThetaBurst')
    ### next set relates to whether a compound EPSP is observed
    '''parser.add_argument('PSPend', type=float, help='expected end of a PSP, typically between 0.265-0.28')
    parser.add_argument('peak2exists', type=bool, help='True if compound/multi EPSP')
    parser.add_argument('pp2exists', type=bool, help='always false')
    parser.add_argument('PP2start', type=float, help='manually define')
    #
    '''
    #parameters only used in GrpAvg
    if flag:
        parser.add_argument('IDfile', type=str)
        parser.add_argument("-samp_time", nargs="+",default=[20,30]) 
        parser.add_argument('-plot_ctrl',type=str,default='100', help='1st bit: show plots, 2nd bit: #columns in figure, 3rd bit: plot correlations') 
        parser.add_argument("-sepvarlist", nargs="+",default=['Status','celltype'],help='list of separation variables for grouping data')
        parser.add_argument("-induction", type=str ,help="induction frequency, e.g. 20 Hz or 10.5 Hz") 

    try:
        args = parser.parse_args(commandline) # maps arguments (commandline) to choices, and checks for validity of choices.
    #if arguments are mapped incorrectly, python wants to exit, but the next line says "don't", instead check whether we are in python (do_exit=False) then don't exit, just give us a warning
    except SystemExit:
        if do_exit:
            raise # raise the exception above (SystemExit) b/c none specified here
        else:
            raise ValueError('invalid ARGS')
    return args