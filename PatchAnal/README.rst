**1. PatchAnal.py**
------------------------

Analyze post-synaptic potentials before and after LTP induction from whole cell patch clamp experiments. Input must be hdf5 format file collected with dPatch.
Assumes that Headstage1 voltage is channel S2, current is channnel S1; headstage2 voltage is channel S4, current is channel S3.  
If a notebook file was saved, it is used to extract routine time and sweep time, relative to breakin, and also read in meta data, such as animal age and ID, tissue region, batah solution.

Analysis includes:

1. Find the PSP in each trace as the most positive peak (after stim time -which is either in notebook file or specified on command line).  

	a. Calculate PSP amplitude relative to resting potential
	
	b. normalize PSP amplitude - divide by mean PSP amplitude determined from pre-induction traces
	
	c. Calculate access resistance from each trace
	
2. IV/IF: from series of current injection, calculate the change in voltage, the number of spikes (if they occur) rheobase, latency to spike, and spike characteristics.  Also calculate rectificiation if trace doesn't have spikes.

3. Extract PSP amplitude from traces used to construct IO curve

4. Analyze the induction protocol.  Determine frequency of within and between bursts, and time between trains.  Also determine whether spikes occur during induction, and the time of spikes relative to the stimulation

Required inputs: experiment name, headstage(s) and celltype(s). If sex, age, drugs, genotype, region are not in the notebook file, these must also be specified. Specify inputs from within python with the following syntax:

ARGS="exper_name -headstage H1 celltype D1-SPN"

or

ARGS="exper_name -headstage H1 celltype D1-SPN -sex F -age 56 -genotype wt -region DL"

*The following additional parameters are specified in the arguments:*

	a. APthresh: minimum amplitude to be considered a spike, default = 0

	b. max_risetime: AP rejected if risetime slower than this, default = 0.004 s

	c. min_risetime: AP rejected if risetime faster than this, default = 0.0003 s

	d. threshval: AP threshold defined as dVm/dt when rise exceeds this fraction of max dV/dt

	e. refract: minimum time between AP, default=0.004 s

	f. basestart: time PRIOR to event for assessing resting membrane potential

	g. base_dur: duration of baseline period, must be smaller than basestart, default=0.10 s

	h. ss_dur: calculate steady state Vm using last X sec of pulse, default=0.050 s

	i. IOrange: values of dig stim for IOtest (calculating IO curve)

	j. induction: name of induction protocol, default = ThetaBurst

	k. inputDir: name of directory containing datafiles for analysis

	l. outputdir: name of directory for output npz files

**2. synaptic.py**
------------------------

Analyze extract and process the electrophysiological recordings of spontaneous PSCs. Input must be hdf5 format file collected with dPatch. Assumes that Headstage1 voltage is channel S2, current is channnel S1.  

Event detection uses the scipy function find_peaks, applied to traces that are lowpass filtered using a 4 pole butterworth filter.  Uses (hard-coded) width and distance parameter value of 0.005 sec in find_peaks for detecting PSCs.

Analysis of detected PSCs includes:

	a. amplitude
	
	b. decay time
	
	c. 10-90% risetime
	
	d. halfwidth (width at half height)
	
	e. inter-event interval
	
Output file contains mean and std of all measures above, experimental parameters from the notebook file, a list of the events, and cumulative density functon of amplitude and IEI.

Required inputs: experiment name, headstage(s) and celltype(s), start and end sweep numbers. Animal ID needs to be in the notebook file. Specify inputs from within python with the following syntax:

ARGS="exper_name -headstage H1 celltype D1-SPN"

*The following additional parameters are specified in the arguments:*

	a. filter: cutoff frequency of butterworth filter, set to 0 to avoid filtering, default=1000Hz

	b. height: prominence parameter in find_peaks for detecting PSCs.  default = 13 pA 
	
	c. wlen:  wlen parameter in find_peaks for detecting PSC, default = 0.1
	
	d. base: duration of time for measuring baseline current - for calculating PSC amplitude, default = 0.004 s
	
	e. decay: duration of time for searching for the 1/e decay point, default=0.03 s

**3. PatchAnalPlots.py**
------------------------

Functions for plotting the following, used in synaptic.py and/or PatchAnal.py

	a. PSP traces - for visual inspection to ensure proper detection and measurement of PSPs
	
	b. Summary of PSP amplitude versus time
	
	c. IV-IF curve traces
	
	d. Induction traces, with artifcats and AP annotated
	
	e. IO curves
	
	f. Histogram of measurements of spontaneous PSCs
	
	g. correlation between two measurements of spontaneous PSCs
	
	
**4. GrpAvgPatchClass.py**
--------------------------

Analyzes groups of experiments - the output of PatchAnal.py. Generates graphs and two types of output data:

    a list of experimental parameters and summary measures, 1 line per experiment, to be used for statistical anaysis
	
    a file of mean, stdev, and N for each group to be used to generate publication quality figures.

*The following parameters are specified in the arguments:* 
	a. IDfile - csv file containing animal information such as sex or genotype (required)

    b. outputdir: full path to location of pickle files (i.e., output files from PatchAnal.py)
	
	c. slope_std_factor: excludes data files in which baseline slope exceeds +/- this factor times the std of the fit to the baseline. default value = 2

    d. sepvarlist: A list of variables and values to used to separate all the data into groups. E.g. ['sex', 'genotype' ] is a list with two separation variables: sex, and genotype.  These variables are used by the Pandas function groupby. The order of specifying variables only matters to how the plots are grouped.
	
	e. samp_time: set of follow-up times (in minutes) for calculating mean plasticity change
	
	f. plot_ctrl: a 3 bit string controlling the plots.

		1. 1st bit: 1 to have only one column of plots, 0 to have one column for each category in second variable of sepvarlist
		
		2. 2nd bit: 1 to plot PSP vs time for each experiment in a group, 0 otherwise
		
		3. 3rd but: 1 to plot correlations between LTP (at summarytime) and age or baseline epsp, 0 otherwise
		
*hard coded parameters:*

	a. nan_threshold: exclude data files if more than this many nans in PSP amplitude. value=10
	
	b. slope_threshold: exclude data files in which baseline slope exceeds +/- this value. value = 0.01 mV/min
	
	c. minimum sweeps: exclude data files if duration of experiment was shorter than this many follow-up sweeps. value = 10

Specify parameters from within python with the following syntax:

    ARGS = "AnimalInfo -outputdir ../pickle/ IDfile -sepvarlist sex region theta -plot_ctrl 000"


**5. GrpAvgSyn.py**
-------------------
Analyzes groups of experiments - the output of synaptic.py. Generates graphs and two types of output data:

    a list of experimental parameters and summary measures, 1 line per experiment, to be used for statistical anaysis
	
    a file of cumulative distribution function (quantiles and cumprob) for each group to be used to generate publication quality figures.

*The following parameters are specified in the arguments:* 
	a. IDfile - csv file containing animal information such as sex or genotype (required)

    b. subdir: full path to location of pickle files (i.e., output files from synaptic.py)
	
    c. sepvarlist: A list of variables and values to used to separate all the data into groups. E.g. ['sex', 'genotype' ] is a list with two separation variables: sex, and genotype.  These variables are used by the Pandas function groupby. The order of specifying variables only matters to how the plots are grouped.
	
	d. plot_ctrl: a 2 bit string controlling the plots.

		1. 1st bit: whether to show graphs or not
		
		2. 2nd bit: 1 to plot correlations between PSC characteristics, 0 otherwise
		
Specify parameters from within python with the following syntax:

    ARGS = "AnimalInfo -subdir ../pickle/ IDfile -sepvarlist sex genotype -plot_ctrl 11"

**6. GrpPlotUtil2.py**
----------------------

Used by GrpAvgPatchClass, to determine mean and stderr versus time for each group, for reading the IDfile, and for group plots

Used by GrpAvgSyn for reading the IDfile

**7. spike_utilities.py**
-------------------------
Functions for extracting spike width, AHP and rectification.  Used in PatchAnal.py

**8. ArgParser.py**
-------------------

argparser used in PatchAnal.py and GrpAvgPatchClass.py
