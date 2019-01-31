from __future__ import print_function as _, division as _
import numpy as np

#in all of these functions, items or times is the spikeTime array (or shuffled array)
def histogram(items, bins):
    """Count items in [[bins[0], bins[1]),
                       [bins[0], bins[1]),
                       ...
                       [bins[0], bins[1])]
    (note the asymmetry)
    """
    bins = bins[:, None]
    fitting = (items >= bins[:-1]) & (items <= bins[1:])
    sum = fitting.sum(axis=1, dtype=histtype)
    assert sum.sum() == (-np.isnan(items)).sum()
    return sum

def histogram2(items, bins):
    left, right = bins[0], bins[-1]
    step = (right - left) / bins.size
    where = (items - left) // step
    assert (where >= 0).all()
    assert (where < bins.size - 1).all()
    hist = np.zeros(bins.size - 1)
    for w in where:
        hist[w] += 1

def histogram3(items, bins):
    left, right = bins[0], bins[-1]
    step = (right - left) / bins.size
    out = np.empty_like(items, dtype=int)
    where = np.floor_divide(items - left, step, out=out)
    k = (where >= 0) & (where < bins.size-1)
    #assert (where >= 0).all()
    #assert (where < bins.size - 1).all()
    return np.bincount(where[k], minlength=bins.size - 1)

def histogram4(items, weight, bins):
    left, right = bins[0], bins[-1]
    step = (right - left) / (bins.size - 1)
    #print('hist4 bins', left,right,step, 'diffs shape', np.shape(items))
    out = np.empty_like(items, dtype=int)
    #the next two steps might be equivalent to k=(items>=left)&(items<right)
    #i.e., only include time differences specified by binMin and binMax
    where = np.floor_divide(items - left, step, out=out, casting='unsafe')
    k = (where >= 0) & (where < bins.size-1)
    #    assert (where >= 0).all()
    #    assert (where < bins.size - 1).all()
    return np.bincount(where[k], weight[k], minlength=bins.size - 1)

def diff_and_hist(args):
    times, i, bins=args
    diffs = times[:, None, i:i+1, :] - times[:, :, None, :]
    diffs[:, i, :, :] = np.nan
    diffs = diffs[-np.isnan(diffs)]
    return histogram3(diffs, bins)

def diff_and_hist2(args):
    times, i, bins=args
    out = np.zeros(bins.size-1)
    for j in range(times.shape[1]):
        if j != i:
            diffs = times[:, None, i:i+1, :] - times[:, j:j+1, None, :]
            diffs = diffs[-np.isnan(diffs)]
            hist = histogram3(diffs, bins)
            norm = hist.sum()
            if norm != 0:
                out += hist / norm
    return out

def diff_and_hist3(args):
    times, i, bins, weights, useFirstInHistogram, useLastInHistogram = args
    #1st index is 5 shuffles, i is the neuron, : is all the spike times
    #finds the time difference between ith neuron and all other neurons for each shuffle
    #using None does a broadcast
    diffs = times[:, None, i:i+1, :, None] - times[:, :, None, None, :]
    # offsets = np.random.rand(*diffs.shape[:3]) * totalTime*2 - totalTime
    # diffs += offsets[..., None, None]
    print('diff&hist',np.shape(times[:,None,i:i+1, :, None]),np.shape(times[:, :, None, None, :]),np.shape(diffs))
    #next, set the difference = np.nan for the ith neuron - eliminate all those spurious zeros
    diffs[:, i, ...] = np.nan
    # diffs[diffs > totalTime] -= 2*totalTime
    # diffs[diffs < -totalTime] += 2*totalTime
    if not useFirstInHistogram: 
        diffs[:, :, :, 0, :] = np.nan
        diffs[:, :, :, :, 0] = np.nan
    if not useLastInHistogram:
        raise NotImplemented
    
    #identify which locations have nans
    k = ~np.isnan(diffs)
    #broadcast the weights to the diffs array
    weight = np.broadcast_arrays(diffs, weights[i][:, None, None, None])[1]
    #extract non-zero diffs and weights in flattened array to compute histogram
    diffs = diffs[k]
    weight = weight[k]
    print('before hist4, i=', i, np.shape(diffs))
    return histogram4(diffs, weight, bins)
