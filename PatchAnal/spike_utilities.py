#utilities for spike properties
import numpy as np

def spike_width(y,x,peaks,height,thresh):
    width=[]
    halfheight=(np.array(height)/2)+thresh
    for i, peakpt in enumerate(peaks):
        beg = end = peakpt
        while beg > 1 and y[beg - 1] > halfheight[i]:
            beg -= 1
        while end + 2 < y.size and y[end + 1] > halfheight[i]:
            end += 1
        width.append( (x[end] + x[end+1] - x[beg] - x[beg-1]) / 2)
    return width

def spike_ahp(y,peaks,widths,inject_endpt,thresh):
    ahp=[]
    ahp_pt=[]
    for i,peakpt in enumerate(peaks):
        end = peaks[i + 1] if i < len(peaks)-1 else inject_endpt
        beg = min(peakpt + int(widths[i]/2), end)
        if end>beg:
            ahp.append( y[beg:end].min()-thresh[i])
            ahp_pt.append(y[beg:end].argmin()+beg)
        else:
            ahp.append(np.nan)
            ahp_pt.append(np.nan)
    return ahp, ahp_pt

def rectification(traces,inject, inj_start,inj_end,window):
    minval=np.zeros(np.shape(inject))
    positions=np.argmin(traces[inj_start:inj_end],axis=0)+inj_start
    for num,pos in enumerate(positions): 
        if inject[num] < 0:
            minval[num]=np.mean(traces[pos-window:pos+window,num])
    return minval

