from __future__ import print_function as _, division as _
import numpy as np
from matplotlib import pyplot as plt
import histograms as hist

def isi_histogram_plot(middles, filename,**hists):
    f = plt.figure()
    f.suptitle('isi histograms, no correction ' + filename)
    f.canvas.set_window_title('isi histogram ' + filename)
    for label, hist in hists.items():
        f.gca().plot(middles, hist, label=label)
    f.gca().set_xlabel('delta T (sec)')
    f.gca().legend()
    f.show()
    return f

def isi_histogram_diff_plot(middles,  filename, hist_diff, ptile95):
    f = plt.figure()
    ptile=ptile95*np.ones(len(middles))
    f.canvas.set_window_title('isi histogram difference (shuffle corrected) ' + filename)
    f.suptitle('isi histogram difference (shuffle corrected) ' + filename)
    ax = f.gca()
    ax.plot(middles, hist_diff, label='unshuffled - shuffled')
    ax.plot(middles, ptile, 'r--', label='95% confidence?')
    ax.set_xlabel('delta T (sec)')
    ax.legend()
    f.show()
    return f

def isi_histogram_by_number_plot(bins, T,  filename, num_neur, weight, additive=False):
    middles = (bins[:-1] + bins[1:]) / 2
    f = plt.figure()
    f.canvas.set_window_title('isi histogram by number ' + filename)
    base = np.zeros_like(middles)
    #put this code (not the plotting part) into function
    for i in range(num_neur):  #loop over k in range(num_neur) for cross-corr
        #Change T[:,i] to T[i,:] to get spikes for single neuron
        spikes = T[:, i][-np.isnan(T[:, i])]
        #spike2=T[:, k][-np.isnan(T[:, k])]
        diffs = spikes[:, None] - spikes #autocorrelations
        #diffs = spikes[:, None] - spikes2 #cross-correlations
        for j in range(diffs.shape[0]):
            diffs[j, j] = np.nan
        histg = hist.histogram4(diffs.ravel(), weight, bins)
        #now, need to shuffle correct each neuron's correlations
        f.gca().plot(middles, base + histg, label='spike #{} ({})'.format(i + 1, j))
        if additive:
            base += histg
    f.gca().legend()
    f.show()
    return f

def isi_histogram_triangle_diff_plot(bins,  filename, histun):
    middles = (bins[:-1] + bins[1:]) / 2
    f = plt.figure()
    f.canvas.set_window_title('isi histogram difference ' + filename)
    ax = f.gca()
    left, right = bins[(histun>0).argmax()], bins[-(histun>0).argmax()]
    total = histun.sum()
    triangle = np.maximum(np.minimum(middles - left, right - middles), 0)
    triangle *= total.sum() / triangle.sum()
    ax.plot(middles, histun - triangle, label='unshuffled - triangle')
    ax.legend()
    f.show()
    return f

def spike_ini_plot(T,  filename, spikers, first=True):

    f = plt.figure()
    if first:
        k = T[:,0].argsort()
    else:
        k = np.empty_like(T[:, 0])
        for i in range(k.size):
            tt = T[i]
            try:
                k[i] = tt[-np.isnan(tt)][-1]
            except IndexError:
                k[i] = np.nan
        k = k.argsort()
    ax = f.gca()
    dat = np.ones(T[0].shape)
    for i in range(T.shape[0]):
        ax.plot(T[k[i]], i*dat, '|')
    f.canvas.set_window_title('spike times ' + filename)
    f.suptitle('spike times ' + filename)
    ax.set_xlabel('time / s')
    ax.set_ylabel('first spike time order')
    ax.set_title('{} spikes'.format(spikers))
    f.show()
    return f

def isi_by_spike_number_plot(T,  filename, maxspike=4):
    f = plt.figure()
    ax = f.gca()
    dT = np.diff(T[:, :maxspike + 1], axis=1)
    for i in range(dT.shape[1]):
        x = dT[:, i]
        ax.hist(x[-np.isnan(x)], alpha=max(0.3, 1-i/30), normed=1, histtype='step', bins=200, label=str(i))
    ax.legend(loc='best')
    f.canvas.set_window_title('isi by spike number (up to {}) {}'.format(maxspike, filename))
    f.show()
    return f

def isi_plot(T,  filename, bins=200):
    f = plt.figure()
    ax = f.gca()
    dT = np.diff(T, axis=1)
    ax.hist(dT[-np.isnan(dT)], normed=1, bins=bins)
    f.canvas.set_window_title('isi ' + filename)
    f.suptitle('isi ' + filename)
    f.show()
    return f

def isi_recurrence_plot(T, filename):
    dT = np.diff(T, axis=1)
    f = plt.figure()
    ax = f.gca()
    ax.scatter(dT[:-1:2], dT[1::2], marker='.', s=1)
    f.canvas.set_window_title('isi recurrence ' + filename)
    f.show()
    return f

def isi_recurrence_plot2(T, filename):
    dT = np.diff(T, axis=1)
    f = plt.figure()
    ax = f.gca()
    x, y = dT[:-1], dT[1:]
    k = ~np.isnan(x) & ~np.isnan(y)
    x, y = x[k], y[k]
    hist, binsx, binsy = np.histogram2d(x, y, bins=800, normed=True)
    im = ax.imshow(np.log(hist), interpolation='none', origin='lower',
                   extent=(binsx.min(), binsx.max(), binsy.min(), binsy.max()))
    f.colorbar(im, shrink=0.5, aspect=10)
    f.canvas.set_window_title('isi recurrence ' + filename)
    f.suptitle('isi recurrence ' + filename)
    f.show()
    return f
