==========
ephys anal
==========
**1. PopSpikeAnal.py**
------------------------
Analyzes individual experiments in which population spikes are recorded before (baseline) and after (follow-up) LTP induction.  The LTP induction is detected automatically as the last gap between pop-spikes greater than induction_gap (e.g. 32 sec, assuming pop-spikes measured once per 30 sec).  Assumes input data is a text file output from labview, in which there are 24 header lines, and the remaining data are in two columns: time and voltage.  Program takes the single voltage column and cuts it into N traces of K samples, where samples_per_trace is read from one of the header lines.

Must provide information on experiments, such as filename that data is in, sex, age, any drugs, brain region, and frequency of induction protocol. Specify these from within python with the following syntax:

   ARGS="filename sex age drugs frequency region"

To modify the information stored with experiments, modify parse_args in pop_spike_utilities.py (and then modify GrpAvgPopSpike.py which uses the same arg parser)

*Adjustable parameters include:*

A. datadir: full path to location of data files
B. outputdir: relative path, relative to python program, for output files,
C. artdecaytime: artifact decay time - if not long enough, then the artifact may be detected as one of the response measures.  If too long, could also miss one of the response measures.
D. artifactthreshold: artifact threshold - artifact is found by looking for the earliest signal above this value.
E. baseline_minutes: how many minutes used to estimate baseline popspike amplitude
F. tracesPerMinute: how often you stimulate and sample per minute
G. sample_times: set of follow-up times (in minutes) for providing mean plasticity change
H. sample_window: number of points either side of sample_times used to calculate mean
I. first_peak_end_fraction: the latest within the trace (units are fraction of trace, from 0 to 1) to look for the main popspike peak
J. slope_std_factor: prints a warning of the baseline slope exceeds +/- this factor times the std of the fit to the baseline.  Make this value large to disable this feature.
K. big_popspike_factor: prints a warning if popspike is greater than this times the mean baseline amplitude.  Also, stores np.nan instead of the amplitude.  Make this value very large (e.g. 10-100) to disable this function

**2. pop_spike_utilities.py**
--------------------------------
Used by PopSpikeAnal.py.  Contains function for plotting the data.  Each trace is plotted along with location of the peak, location of positive peak, location where search for the pop-spike and fiber volley begin.  This allows user to verify the data extraction, and determine if parameters such as artifact decay time need to be adjusted.

**3. PSPanalSA.py**
------------------------
Analyze post-synaptic potentials before and after LTP induction from whole cell patch clamp experiments. Assumes input data are a set of files in igor binary.

Must provide information on experiments, such as filename that data is in, sex, age, any drugs, celltype, genotype, light level, light responsive (the latter two are for optogenetic protocols), depolarization (either via somatic current injection or optogenetically), whether APs were observed during the induction protocol, time between break-in and beginning of LTP induction, tip resistance, temperature of bath.  Specify these from within python with the following syntax:

   ARGS="20-Aug-2015_SLH001 M 29 heat nodrug 5.4 12 A2a+ MSN non 0 soma APs"

To modify the information stored with experiments, modify parse_args 

*Additional parameters (common to all experiments):*

a. hyperstart,hyperend, Iaccess: start time, end time and current amplitude of hyperpolarizing pulses injected to monitor series resistance
b. basestarttime, baseendtime: time period prior to electrical stimulation to use for membrane potential for calculating PSP amplitude
c. dt: interval between samples (1/sampling frequency)
d. inputDir: name of directory containing datafiles for analysis.

*File locations:*

a. FileDir=inputDir+args.experiment+"_Waves/" specifies location of data files for an individual experiment
b. outputDir: name of directory to write output files
c. filenameending: text string appended to name of experiment which specifies wildcard pattern filenames: pattern=FileDir+"W"+args.experiment+filenameending

**4. GrpAvgPopSpike.py**
------------------------
Analyzes groups of experiments - the output of PopSpikeAnal.py.
Generates graphs and two types of output data:

A. a list of experimental parameters and summary measures, 1 line per experiment, to be used for statistical anaysis
   
B. a file of mean, stdev, and N for each group to be used to generate publication quality figures.
   
*Adjustable parameters include:*

1. subdir: full path to location of pickle files (i.e., output files from PopSpikeAnal.py)
2. slope_std_factor: Currently not used.  Could be used to excludes data files in which baseline slope exceeds +/- this factor times the std of the fit to the baseline.  Instead, we are using ...
3. slope_threshold: exclude data files in which baseline slope exceeds +/- this value.
4. minimum_sweeps: Excludes data files which have insufficient traces following the induction.  This number is total number of minutes, which user calculates from baseline_minutes + follow-up minutes
5. sample_times: set of follow-up times (in minutes) for providing mean plasticity change
6. sepvarlist: A list of variables and values to used to separate all the data into groups. E.g.[ ['sex',['F','Fe','M']], ['drug', ['none']] ] is a list with two separation variables: sex, which can have one of 3 values, and drug, which could have multiple values, but by indicating a single value, the code will give two drug groups: none, and everything else. The order of specifying variables only matters to how the plots are grouped.  

**5. GrpAvgPatchMultiGroups.py**
---------------------------------
Analyzes groups of experiments - Similar to GrpAvgPopSpike, but uses output of PSPanalSA.py.
Generates graphs and two types of output data:

A. a list of experimental parameters and summary measures, 1 line per experiment, to be used for statistical anaysis
  
B. a file of mean, stdev, and N for each group to be used to generate publication quality figures.
   
*Adjustable parameters include:*

A. subdir: full path to location of pickle files (i.e., output files from PSPanalSA.py)
B. sepvarlist: A list of variables and values to used to separate all the data into groups - see explanation under GrpAvgPopSpike.py
C. meanstart,meanend: trace numbers corresponding to 15-20 min after induction

This program was used in Hawes et al. J Physiology 2017 for analysis of data.

**6. GrpPlotUtil.py**
-------------------------
Used by GrpAvgPopSpike.py and by GrpAvgPatchMultiGroups.py 

**7. TBSanal.py**
-------------------------

**8. AnalyzeIV.py**
-------------------------
Analyzes IF and IV curves from whole cell patch clamp experiments.
Assumes IF is separate set of curves from IV.  Must specify (or use default values) or starting current injection and increment.  Must specify (or use default values) for time of current injection onset and duration of current injection.

**9. HVAanal.py**
-------------------------
Analyze two pulse voltage clamp experiments from whole cell patch clamp experiments to determine calcium dependent inactivation of calcium currents.

Assumes input data is are a set of files in igor binary.

This program was used in Evans et al. J Neurophysiology 2015 for analysis of data.

**10. RampAnal.py**
-------------------------
Analyze ramp voltage clamp from whole cell patch clamp experiments in order to extract leak conductance.

Assumes input data is are a set of files in igor binary.

This program was used in Evans et al. J Neurophysiology 2015 for analysis of data.

**11. SASdataIF.py**
-------------------------

**These python programs contain utilities used by TBSanal.py and AnalyzeIV.py for spike dection and characterization:**
1. compat.py
2. detect.py
3. loader.py
4. signal_smooth.py
5. utilities.py
