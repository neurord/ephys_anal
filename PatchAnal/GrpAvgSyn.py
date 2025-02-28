import numpy as np
import sys  
import pandas as pd 
import glob
import os
print (os.getcwd())
from matplotlib import pyplot
import argparse as argp
from GrpPlotUtil2 import group_to_word

def ArgParserSyn(commandline,do_exit):
    parser = argp.ArgumentParser()
    parser.add_argument('-datadir',type = str, default = "C:\\Users\\klblackwell\\Documents\\Python\\ephys_anal\\PatchAnal\\IPSCs\\",   help="Directory containing raw data")
    parser.add_argument('-plot_ctrl',type=str,default='11', help='1st bit: show plots, 2nd bit: plot correlations') 
    parser.add_argument("-sepvarlist", nargs="+",default=['sex','region'],help='list of separation variables for grouping data')
    #specify these next variables to analyze only a subset of the data
    parser.add_argument("-drug", type=str, choices=["ctrl","CNQX","Picro"],help="ctrl, CNQX or Picro")
    parser.add_argument("-region", type=str, choices=["DM", "DL"], help="which striatal region") 
    parser.add_argument("-celltype", type=str, choices=["D1", "D2"], help="which celltype") 
    try:
        args = parser.parse_args(commandline) # maps arguments (commandline) to choices, and checks for validity of choices.
    #if arguments are mapped incorrectly, python wants to exit, but the next line says "don't", instead check whether we are in python (do_exit=False) then don't exit, just give us a warning
    except SystemExit:
        if do_exit:
            raise # raise the exception above (SystemExit) b/c none specified here
        else:
            raise ValueError('invalid ARGS')
    return args

class GrpSyn:
    def __init__(self, params): 
        self.params = params
        self.datadir = params.datadir
        self.sepvarlist= params.sepvarlist
        self.plot_corr=int(params.plot_ctrl[1]) #to plot correlation between LTP and age or baseline epsp

    def read_data(self): 
        DATAS = []
        PARAMS = []
        for outfname in self.outfnames:
            with open(outfname, 'rb') as f:
                datadict = np.load(f,allow_pickle=True)
                data=datadict['data'].item()
                param=datadict['params'].item()
            ignore = ((self.params.drug and self.params.drug != param['drug']) or
                        (self.params.region and self.params.region != param['region']))
                    # or(self.params.sex and self.params.sex != param['sex'])                        
            if ignore: #regardless of whether data meets exclusion criteria
                if self.print_info:
                    print ("***** ignoring:", param['exper'])
                next 
            else: #data meets all criteria
                PARAMS.append(param)
                DATAS.append(data)                    
        dfdata = pd.DataFrame(DATAS)
        dfparams=pd.DataFrame(PARAMS)
        self.whole_df=pd.concat([dfparams,dfdata], axis = 1)

    def group_data(self): 
        self.grp_data = self.whole_df.groupby(self.sepvarlist)
        self.df_num=self.whole_df.select_dtypes(exclude=[object])
        agg={}
        for col in self.df_num.columns:
            agg[col]='mean'
        self.mean_df=self.grp_data.agg(agg)

    def plot_groups(self,ycol):
        #grp_nm=group_to_word(group)
        xvals=self.grp_data.groups
        x=[group_to_word(g) for g in xvals.keys()]
        yvals=self.grp_data.mean(numeric_only=True)[ycol]
        pyplot.scatter(x,yvals)
        pyplot.ylabel(ycol)

    def write_stat_data(self):
        self.df_str=self.whole_df.select_dtypes(include=[object])
        params = self.df_str.values #into a 2d array you write it into SASoutput.
        values=[]
        for ii in range(len(self.whole_df)):
            row=self.whole_df.iloc[ii][self.df_num.columns].values
            newrow= [str(np.round(val,6)) for val in row]
            values.append(newrow)
        SASoutput=np.hstack((params,np.array(values)))
        SASheader= '   '.join(self.whole_df.columns)
        f=open(self.datadir+"PARAMSforSTATS.txt", 'w')
        f.write(SASheader +"\n")
        np.savetxt(f, SASoutput, fmt='%s', delimiter='   ')
        f.close()


if __name__ =='__main__':        
    ARGS = "-sepvarlist region -plot_ctrl 11"      
    try:
        commandline = ARGS.split() 
        do_exit = False
    except NameError: 
        commandline = sys.argv[1:]
        do_exit = True
    params=ArgParserSyn(commandline,do_exit)
    print(commandline, params)
    grp=GrpSyn(params)
    grp.pattern = grp.datadir+'*.npz'
    grp.outfnames = sorted(glob.glob(grp.pattern))
    if not len(grp.outfnames):
        sys.exit('************ No files found *********')
    if len(grp.outfnames):
        grp.read_data()
        grp.group_data()
        grp.write_stat_data()
        grp.plot_groups('decay_mean')
 