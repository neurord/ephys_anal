import numpy as np
import glob
import pickle
import csv

######################## How to add a new parameter to pickle files #########
directory='C:/Users/vlewitus/Documents/Python Scripts/Pickle/'
pattern = directory+'*.pickle'
outfnames = glob.glob(pattern)
E2_value_file=directory+'E2values.csv'
file_E2_dict={}
with open(E2_value_file) as f: #only needed to open file
	csvreader = csv.reader(f, delimiter=',')
	next(csvreader) #skips the header
	for row in csvreader:
		file_E2_dict[row[0]]=row[1]
        #print(row) ['f1','3.56']
#this has to be run once to give estradiol variable to every file
#will work as long as E2_value_file has every estradiol measurement made
for outfname in outfnames:
	with open(outfname,'rb') as f:
		datadict = np.load(f,fix_imports=True)
		if datadict['parameters'].exper in file_E2_dict.keys():
			datadict['parameters'].estradiol=float(file_E2_dict[datadict['parameters'].exper])
			print(datadict['parameters'])
		#elif datadict['parameters'].estradiol >-1
			#datadict['parameters'].estradiol = -1
	with open(outfname,'wb') as f:
		pickle.dump(datadict,f)

'''
#continued updates
#will work with incremental files of new values or file with all values, but won't add estradiol to every other file
for fname in file_E2_dict.keys():
	outfname=directory+fname+'.pickle'
	with open(outfname,'rb') as f:
		datadict=np.load(f,fix_imports=True)
		datadict['parameters'].estradiol=file_E2_dict[fname]
	with open(outfname,'wb') as f:
		pickle.dump(datadict,f)
'''