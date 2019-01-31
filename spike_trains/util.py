from __future__ import print_function as _, division as _
import numpy as np
import csv

def readcsv(fname):
    data=[]
    with open(fname) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            data.append(row)
    neurnames=[name for name in data[0] if len(name)>0]
    numrows=len(data[1:])
    numcols=len(neurnames)
    print('file,num_neurons',fname, numcols, 'names',neurnames)
    timestamp=np.empty((numrows,numcols))
    for index,row in enumerate(data[1:]):
        timestamp[index]=[float(row[j]) if len(row[j])>0 else np.nan for j in range(numcols) ]
    return timestamp

    def trains():
        meantrains = 3
        meanspikes = 1.5
        intertrainisi = 0.012
        intertrainisinoise = 0
        intertrainnoise = 0.001
        min_isi = 0.004

        time = 0
        while time < binMax:
            t = np.random.exponential(binMax / meantrains)
            time += t

            n = np.random.geometric(1 / meanspikes)
            itisi = intertrainisi + np.random.randn(n) * intertrainisinoise
            isis = np.cumsum(min_isi + np.random.exponential(itisi - min_isi, size=n))
            times = time + isis
            yield times[(times < binMax) & (times >= 0)]
            time = times[-1]

