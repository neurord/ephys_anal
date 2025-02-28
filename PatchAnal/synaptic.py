import numpy as np
import sys
import argparse
import pandas as pd
from spike_utilities import find_notebook_file
from PatchAnalPlots import rows_columns

HEADSTAGE_I={'S1': 'H1', 'S3': 'H2'}

def ArgParserSyn(commandline,do_exit):
    parser = argparse.ArgumentParser()
    parser.add_argument('exper', type=str, help = 'give name of experiment')
    parser.add_argument('-datadir',type = str, default = "C:\\Users\\klblackwell\\Documents\\Python\\ephys_anal\\PatchAnal\\IPSCs\\",   help="Directory containing raw data")
    parser.add_argument('-headstages', type=str, nargs="+", help='which headstage to analyze (default=both)', default=['H1','H2']) #FIXME - make this required for PatchAna?
    parser.add_argument("-celltypes", type=str, nargs="+", choices=["D1-SPN", "FSI", "D2-SPN",'ChI','other'],help="D1-SPN, D2-SPN, FSI or other cell type, one for each headstage") #FIXME - make this required for PatchAnal?
    try:
        args = parser.parse_args(commandline) # maps arguments (commandline) to choices, and checks for validity of choices.
    #if arguments are mapped incorrectly, python wants to exit, but the next line says "don't", instead check whether we are in python (do_exit=False) then don't exit, just give us a warning
    except SystemExit:
        if do_exit:
            raise # raise the exception above (SystemExit) b/c none specified here
        else:
            raise ValueError('invalid ARGS')
    return args

class SynAnal():
    def __init__(self,params): #change to PopSpike(args)
        self.datadir = params.datadir
        #self.outputdir = params.outputdir
        self.exper=params.exper
        if len(params.headstages)==len(params.celltypes): 
            self.celltypes={h:params.celltypes[i] for i,h in enumerate(params.headstages)}
        else:
            print('*********** Need to specify celltype for each valid headstage ***********')
            exit()
        self.params={'exper':self.exper}
        
    def read_notebook(self): #extract routine and sweep time from notebook file, if it was saved
        def parse_lines(all_lines):
            for i,line in enumerate(all_lines):
                if line.startswith('Mini Results'): #start of synaptic data for one headstage
                    lin=line.strip('\n')
                    event_start={lin:i}                    
                    while event_start[lin]<len(all_lines): #keep reading until the first numeric line
                        entry=all_lines[event_start[lin]].split(',')[0].strip(' ')
                        if entry.isnumeric():
                            event_end={lin:event_start[lin]}
                            break
                        else:
                            event_start[lin]+=1
                    for j,ln in enumerate(all_lines[event_start[lin]:]): #Now, keep reading till the end of the data for that headstage
                        if ln.startswith('\n'):
                            event_end[lin]+=j
                            break
                    channel=line.split('_')[-1].strip('\n')
                    headstage=HEADSTAGE_I[channel]
                    celltype=self.celltypes[headstage]
                elif line.startswith('"Entry"'):
                    full_line=''.join([all_lines[j].strip('\n') for j in range(i,event_start[lin])])
                    columns=full_line.split(',')
                    column_names=[c.strip(' ').strip('"') for c in columns]
                elif line.startswith('1 - '): #remaining lines: read meta data
                    if 'Animal Age' in line:
                        self.age=self.params['age']=line.split(', ')[-1][0:-1]
                    elif 'Animal Identifier' in line:
                        self.ID=self.params['ID']=line.split(', ')[-1][0:-1]
                    elif 'Tissue Preparation' in line:
                        self.region=self.params['region']=line.split(', ')[-1][0:-1] #override default
                    elif 'Animal Preparation' in line:
                        self.SurgeryDate=self.params['SxDate']=line.split(', ')[-1][0:-1] #override default
                    elif 'Bath Solution' in line:
                        self.drug=self.params['drug']=line.split(', ')[-1][0:-1] #override default
                #else:
                    #do nothing with remaining lines
            return event_start,event_end, column_names

        import glob
        import os
        nfname=self.datadir+'Notebook_20'+'*.txt'  #find all Notebook files with correct date
        files=glob.glob(nfname)
        if len(files):
            findex,file_found=find_notebook_file(files,self.exper)
            file_found=True
        else:
            file_found=False
        if file_found: 
            print('Using notebook file:',files[findex])
            self.params['file']=os.path.basename(files[findex]).split('.')[0]  
            with open(files[findex],'r') as myfile: 
                all_lines=myfile.readlines()
            event_start,event_end,column_names=parse_lines(all_lines)
            self.dfset={}
            for key in event_start.keys(): #will need event_end when there are two headstages
                PSC_list=[]
                for line in range(event_start[key],event_end[key]):
                    data=all_lines[line].split(',')
                    temp_dict={}
                    for col,entry in zip(column_names,data):
                        temp_dict[col]=float(entry.strip().split()[0])
                    PSC_list.append(temp_dict)
                self.dfset[key] = pd.DataFrame.from_dict(PSC_list) 
        else:
            print('0 or multiple Notebook files found using pattern:',nfname)
        return

    def analyze(self):
        self.results={k:{} for k in self.dfset.keys()}
        self.measures={'IEI':'Interevent Interval (s)', 
                  'rise':'10-90% Rise Time (s)', 
                  'decay':'Event Decay Tau (s)',
                  'amp':'Event Amplitude (A)',
                  'auc':'Event Integral (A*s)'}
        for key in self.dfset.keys():
            for m,nm in self.measures.items():
                self.results[key][m+'_mean']=self.dfset[key].mean()[nm]
                self.results[key][m+'_stdev']=self.dfset[key].std()[nm]
            self.results[key]['Freq_mean']=np.nanmean(1/self.dfset[key]['Interevent Interval (s)'])
            self.results[key]['Freq_std']=np.nanstd(1/self.dfset[key]['Interevent Interval (s)'])

    def plot_hist(self,measures):
        from matplotlib import pyplot as plt
        plt.ion()
        n=len(measures)
        cols,rows=rows_columns(n)
        fig,axes=plt.subplots(rows,cols)
        axes=fig.axes
        for key in self.dfset.keys():
            for ax,m in enumerate(measures):
                minv=np.percentile(self.dfset[key][m].dropna().values,1)
                maxv=np.percentile(self.dfset[key][m].dropna().values,99)
                axes[ax].hist(self.dfset[key][m].values,range=(minv,maxv),bins=20)
                axes[ax].set_xlabel(m)
        fig.canvas.manager.set_window_title('Histogram '+self.params['exper'])
        fig.tight_layout()
        return fig
    
    def plot_pairs(self,pairs):
        from matplotlib import pyplot as plt
        plt.ion()
        cols,rows=rows_columns(len(pairs))
        fig,axes=plt.subplots(rows,cols)
        axes=fig.axes
        for ax,p in enumerate(pairs):
            for key in self.dfset.keys():
                corr=self.dfset[key][p[0]].corr(self.dfset[key][p[1]])
                axes[ax].scatter(self.dfset[key][p[0]],self.dfset[key][p[1]],label=key+' '+str(round(corr,4)))
                axes[ax].set_xlabel(p[0])
                axes[ax].set_ylabel(p[1])
                axes[ax].legend()            
        
    def save_data(self):
        self.params['SxDate']=self.SurgeryDate
        self.params['age']=self.age
        self.params['ID']=self.ID
        self.params['region']=self.region
        self.params['drug']=self.drug
        if 'exper' in self.params:
            outfname=self.datadir+self.params['exper']
        else:
            outfname=self.datadir+self.params['file']
        for key in self.dfset.keys():
            h=HEADSTAGE_I[key.split('_')[-1]]
            self.params['celltype']=self.celltypes[h]
            np.savez(outfname+h, params=self.params,data=self.results[key])


if __name__=='__main__':
    ARGS='250220_4 -headstages H2 -celltype D1-SPN'
    try:
        commandline = ARGS.split() 
        do_exit = False
    except NameError: 
        commandline = sys.argv[1:]
        do_exit = True

    params=ArgParserSyn(commandline,do_exit)
    print(commandline, params)
    exp=SynAnal(params)
    exp.read_notebook()
    for key,df in exp.dfset.items():
        print(key,'           correlations\n',df.corr())
    print()
    exp.analyze()
    exp.save_data()
    fighist=exp.plot_hist(list(exp.measures.values()))
    pairs=[('Event Amplitude (A)','Event Decay Tau (s)'),('Event Decay Tau (s)','10-90% Rise Time (s)')]
    fig=exp.plot_pairs(pairs)
#TODO:
#   - GroupAVG



    
