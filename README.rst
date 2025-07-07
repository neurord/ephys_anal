==========
ephys anal
==========
**1. PopSpikeAnalClass.py**
------------------------
Analyzes individual experiments in which population spikes are recorded before (baseline) and after (follow-up) LTP induction.  The LTP induction is detected automatically as the last gap between pop-spikes greater than induction_gap. Induction gap is a pause 10% greater than time between traces. Assumes input data is a text file output from labview, in which there are 24 header lines, and the remaining data are in two columns: time and voltage.  Program takes the single voltage column and cuts it into N traces of K samples, where samples per trace is read from one of the header lines.

Finds the PopSpike in each trace as the most negative peak after the artifact decay + fiber volley.  Then calculates amplitude of PopSpike relative to either the baseline prior to artifact or the positive peak - most positive value between artifact decay and negative peak.
Calculates the baseline PopSpike as the mean of PopSpikes during the baseline period prior to induction.  Then calculates a normalized PopSpike as ratio of PopSpike amplitude to baseline PopSpike.
Calculates the slope of the baseline PopSpike, and prints warning if the slope is statistically significant.
Graphs all traces with location of PopSpike, baseline and positive peak annotated for visual inspection.
Creates output file with PopSpike, analysis parameters and experimental parameters, for use in GrpAvgPopSpikeClass.

Must provide information on experiments, such as filename that data is in, sex, age, any drugs, brain region, and frequency of induction protocol. Specify these from within python with the following syntax:

   ARGS="filename sex age drugs frequency region"

To modify the information stored with experiments, modify parse_args in pop_spike_utilities.py (and then modify GrpAvgPopSpikeClass.py which uses the same arg parser)

*Adjustable parameters include:*

A. datadir: full path to location of data files
B. outputdir: relative path, relative to python program, for output files,
C. decay: artifact decay time - if not long enough, then the artifact may be detected as one of the response measures.  If too long, could also miss one of the response measures.
D. artthresh: artifact threshold - artifact is found by looking for the earliest signal above this value.
E. base_min: how many minutes used to estimate baseline popspike amplitude
F. samp_time: set of follow-up times (in minutes) for providing mean plasticity change
G. FVwidth: fiber volley width.  Begin search for pop spike after artifact decay + FVwidth
H. samp_win: number of points either side of samp_time used to calculate mean
I. first_peak_end_frac: the latest within the trace (units are fraction of trace, from 0 to 1) to look for the main popspike peak
J. slope_std: prints a warning of the baseline slope exceeds +/- this factor times the std of the fit to the baseline.  Make this value large to disable this feature.
K. angle: angle the slices were cut
L. base_min: number of minutes prior to induction gap for measuring baseline

This program was used in Lewitus and Blackwell ENeuro 2023 for analysis of data.

**2. pop_spike_utilities.py**
--------------------------------
Used by PopSpikeAnal.py.  Contains arg parser and function to read labview file format.

Function for plotting the data.  Each trace is plotted along with location of the peak, location of positive peak, location where search for the pop-spike and fiber volley begin.  This allows user to verify the data extraction, and determine if parameters such as artifact decay time need to be adjusted.

**3. PopSpikeIO.py**
------------------------
Similar to PopSpikeAnalClass, but analyzes the first 30 traces which measures PopSpike versus stimulation amplitude.  Uses analysis parameters from the output of PopSpikeAnalClass.

*Adjustable parameters include:*

A. datadir: full path to location of data files
B. outputdir: relative path, relative to python program, for output files from PopSpikeAnalClass
C. Smallest current and Largest current used for IO curve - assumes values are equally spaced.

This program was used in Lewitus and Blackwell ENeuro 2023 for analysis of data.

**4. GrpAvgPopSpikeClass.py**
------------------------
Analyzes groups of experiments - the output of PopSpikeAnalClass.py.
Generates graphs and two types of output data:

A. a list of experimental parameters and summary measures, 1 line per experiment, to be used for statistical anaysis
   
B. a file of mean, stdev, and N for each group to be used to generate publication quality figures.
   
*Adjustable parameters include:*

1. subdir: full path to location of pickle files (i.e., output files from PopSpikeAnal.py)

2. slope_std_factor: excludes data files in which baseline slope exceeds +/- this factor times the std of the fit to the baseline. default value = 2 

4. minimum_sweeps: Excludes data files which have insufficient traces following the induction.  This number is total number of minutes, which user calculates from baseline_minutes + follow-up minutes

5. sample_times: set of follow-up times (in minutes) for providing mean plasticity change. defaults = [30,60,90,120]

6. sepvarlist: A list of variables and values to used to separate all the data into groups. E.g.[ ['sex',['F','Fe','M']], ['drug', ['none']] ] is a list with two separation variables: sex, which can have one of 3 values, and drug, which could have multiple values, but by indicating a single value, the code will give two drug groups: none, and everything else. The order of specifying variables only matters to how the plots are grouped.  

Specify these from within python with the following syntax:

   ARGS = "-outputdir ../pickle/ -sepvarlist sex region theta -plot_ctrl 000" 

*hard coded parameters:*

1. slope_threshold: exclude data files in which baseline slope exceeds +/- this value. value = 0.01 mV/min
2. nan_threshold: exclude data files if more than this many nans in PopSpike amplitude.  value=20 
3. minimum sweeps: exclude data files if duration of experiment was shorter than this many follow-up sweeps

This program was used in Lewitus and Blackwell ENeuro 2023 for analysis of data.

**5. GrpAvgPopSpikeIO.py**
---------------------------------
Analyzes IO curves for groups of experiments - Similar to GrpAvgPopSpikeClass.
Generates graphs and two types of output data:

A. a list of experimental parameters and summary measures, 1 line per experiment, to be used for statistical anaysis
  
B. a file of mean, stdev, and N for each group to be used to generate publication quality figures.
   
*Adjustable parameters include:*

A. subdir: full path to location of pickle files (i.e., output files from PopSpikeIO.py)
B. sepvarlist: A list of variables and values to used to separate all the data into groups - see explanation under GrpAvgPopSpikeClass.py
C. plot_ctrl: a 3 bit string controlling the plots.

   1. 1st bit: 1 to have only one column of plots, 0 to have one column for each category in second variable of sepvarlist 
   2. 2nd bit: 1 to plot PopSpike vs time for each experiment in a group, 0 otherwise
   3. 3rd but: 1 to plot correlations between LTP (at summarytime) and age or baseline epsp, 0 otherwise

This program was used in Lewitus and Blackwell ENeuro 2023 for analysis of data.

**6. GrpPlotUtil2.py**
-------------------------
Used by GrpAvgPopSpikeClass.py 

**7. addParam.py**
-------------------------

Read in all the PopSpikeAnalClass datafiles, and add one additional experimental parameter from csv file.


