# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 15:32:19 2018

@author: kblackw1
"""
import numpy as np
from matplotlib import pyplot as plt

plt.ion()
f = plt.figure()
ax = f.gca()
mouse='Mouse2'
middles=np.loadtxt(mouse+'_bins.txt')
f.suptitle(mouse)
for condition in ['saline','ethanol']:
    fname=mouse+'_'+condition+'.csv_hist.txt'
    data=np.loadtxt(fname)
    ax.plot(middles,data,label=condition)
ax.set_xlabel('time / s')
ax.legend()