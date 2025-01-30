# -*- coding:utf-8 -*-
from __future__ import print_function, division
from future_builtins import zip
try:
    range = xrange
except NameError:
    pass

import glob
import numbers
import contextlib
import functools
import os
import math
import operator
import copy
from collections import namedtuple
import numpy as np
from scipy import optimize
from igor import binarywave

from compat import TemporaryDirectory

import utilities
from signal_smooth import smooth
from detect import detect_peaks

class vartype(object):
    def __init__(self, x, dev=0):
        self.x = x
        self.dev = dev

    @property
    def positive(self):
        return self.x > self.dev*3

    @property
    def negative(self):
        return self.x < -self.dev*3

    def __nonzero__(self):
        return abs(self.x) > self.dev*3

    def __sub__(self, other):
        return vartype(self.x - other.x, (self.dev**2 + other.dev**2)**0.5)

    def __add__(self, other):
        return vartype(self.x + other.x, (self.dev**2 + other.dev**2)**0.5)

    def __rsub__(self, other):
        return vartype(other - self.x, self.dev)

    def __lt__(self, other):
        return self.x < other.x

    def __truediv__(self, other):
        if isinstance(other, numbers.Real):
            return vartype(self.x/other, self.dev/other)
        else:
            return NotImplemented

    def __str__(self):
        preca = int(math.floor(math.log10(abs(self.x))))
        precb = int(math.floor(math.log10(self.dev)))
        prec = -min(preca, precb, 0)
        return '{0.x:.{1}f}±{0.dev:.{1}f}'.format(self, prec)

    def __repr__(self):
        preca = int(math.floor(math.log10(abs(self.x))))
        precb = int(math.floor(math.log10(self.dev)))
        prec = -min(preca, precb, 0) + 1
        return '{0.__class__.__name__}({0.x:.{1}f}, {0.dev:.{1}f})'.format(self, prec)

    @classmethod
    def average(cls, vect):
        sq = vect.dev**-2
        var = 1 / sq.sum()
        x = (vect.x * sq * var).sum()
        return cls(x, var**0.5)

    @classmethod
    def array(cls, items):
        x = np.array([p.x for p in items])
        dev = np.array([p.dev for p in items])
        return np.rec.fromarrays((x, dev), names='x,dev')

def load_dir(dir, timestep=1e-4):
    files = sorted(glob.glob(os.path.join(dir, '*.ibw')))
    inputs = (binarywave.load(file)['wave']['wData'] for file in files)
    datas = (np.rec.fromarrays((np.arange(input.size)*timestep, input), names='x,y')
             for input in inputs)
    return np.vstack(datas).view(np.recarray)

def array_mean(data):
    return vartype(data.mean(), data.var(ddof=1)**0.5)

def array_diff(wave, n=1):
    xy = (wave.x[:n] + wave.x[n:])/n, np.diff(wave.y)
    return np.rec.fromarrays(xy, names='x,y')

def array_sub(reca, recb):
    xy = (reca.x - recb.x, (reca.dev**2 + recb.dev**2)**0.5)
    return np.rec.fromarrays(xy, names='x,dev')

def array_rms(rec):
    return ((rec.x / rec.dev)**2).mean()**0.5

def _find_baseline(wave, before=.2, after=0.75):
#Returns mean value of wave after excluding "outliers", values > 95th or < 5th percentile
    what = wave.y[(wave.x < before) | (wave.x > after)]
    cutoffa, cutoffb = np.percentile(what, (5, 95))
    cut = what[(what > cutoffa) & (what < cutoffb)]
    return array_mean(cut)

def _find_spikes(wave, min_height=0.0, maximum_spike_charge_time=.004,threshval=0.02):
    peaks = detect_peaks(wave.y)
    peaks = peaks[wave.y[peaks] > min_height]

    thresholds = np.empty(peaks.size)
    for i in range(len(peaks)):
        start = (wave.x >= wave.x[peaks[i]] - maximum_spike_charge_time).argmax()
        x = wave.x[start:peaks[i] + 1]
        y = wave.y[start:peaks[i] + 1]
        yderiv = np.diff(y)
        #spike threshold is point where derivative is 2% of steepest
        ythresh = threshval * yderiv.max()
        thresh = y[yderiv > ythresh].min()
        thresholds[i] = thresh
    return peaks, thresholds

def _find_steady_state(wave, after=.25, before=0.6, cutoff_percentile=80):
#Returns mean value of wave for t between "after" and "before" times
#after excluding "outliers", values > 80th percentile (which is parameter)
#80th percentile excludes the spikes
    data = wave.y[(wave.x > after) & (wave.x < before)]
    cutoff = np.percentile(data, cutoff_percentile)
    cut = data[data < cutoff]
    return array_mean(cut)

def _find_falling_curve(wave, window=20, after=0.2, before=0.6):
    d = array_diff(wave)
    dd = smooth(d.y, window='hanning')[(d.x > after) & (d.x < before)]
    start = end = dd.argmin() + (d.x <= after).sum()
    while start > 0 and wave[start - 1].y > wave[start].y and wave[start].x > after:
        start -= 1
    sm = smooth(wave.y, window='hanning')
    smallest = sm[end]
    # find minimum
    while (end+window < wave.size and wave[end+window].x < before
           and sm[end:end + window].min() < smallest):
        smallest = sm[end]
        end += window // 2
    ccut = wave[start + 1 : end]
    return ccut

simple_exp = lambda x, amp, tau: amp * np.exp(-(x-x[0]) / tau)
negative_exp = lambda x, amp, tau: amp * (1-np.exp(-(x-x[0]) / tau))
falling_param = namedtuple('falling_param', 'amp tau')
function_fit = namedtuple('function_fit', 'function params')

def _fit_falling_curve(ccut, baseline, steady):
    if ccut.size < 5:
        func = None
        params = falling_param(vartype(np.nan, np.nan),
                               vartype(np.nan, np.nan))
    else:
        init = (ccut.y.min()-baseline.x, ccut.x.ptp())
        func = negative_exp if (steady-baseline).negative else simple_exp
        popt, pcov = optimize.curve_fit(func, ccut.x, ccut.y-baseline.x, (-1,1))
        pcov = np.zeros((2,2)) + pcov
        params = falling_param(vartype(popt[0], pcov[0,0]**0.5),
                               vartype(popt[1], pcov[1,1]**0.5))
    return function_fit(func, params)

def _find_rectification(ccut, steady, window_len=11):
    if ccut.size < window_len + 1:
        return vartype(np.nan)
    pos = ccut.y.argmin()
    end = max(pos + window_len//2, ccut.size-1)
    bottom = array_mean(ccut[end-window_len : end+window_len+1].y)
    return steady - bottom

class Params(object):
#These parameters may need to be adjusted, especially steady before and after, and baseline
#Better to change these to parameters determined in the main program
#These are times, in sec used for measurements
#Baseline: calculated from part of trace before "baseline_before" and after "baseline_after
    baseline_before = 0.2
    #current injection stops at 0.6, wait a bit for it to return to steady state
    baseline_after = 0.75

    #(steady_after - baseline_before) indicates how much time to wait after initiating current injection before calculating steady state
    steady_after = .25
    #current injection stops at 0.6
    steady_before = .6
    steady_cutoff = 80

    falling_curve_window = 20
    rectification_window = 11

Fileinfo = namedtuple('fileinfo', 'experiment sarah protocol series sweep')

def _calculate_current(fileinfo, IV, IF,IVseries):
#calculatese the set of current injection values for either IV or IF curve
    print(fileinfo)
    assert fileinfo.protocol == 1
    start, inc = IV if fileinfo.series == IVseries else IF
    return start + inc * (fileinfo.sweep - 1)

class IVCurve(object):
    def __init__(self, filename, fileinfo, injection, x, y, params):
        self.filename = filename
        self.fileinfo = fileinfo
        self.injection = injection
        self.params = params
        self.wave = np.rec.fromarrays((x, y), names='x,y')

    @classmethod
    def load(cls, dirname, filename, IV, IF, params,IVseries, time):
        path = os.path.join(dirname, filename)
        wave = binarywave.load(path)['wave']
        data = wave['wData']
        dtfile = wave['wave_header']['hsA']
        time = np.arange(len(data)) * dtfile

        #what does basename do: it just returns the filename
        #why is this fileinfo important - allows us to refer to aspects of the filename
        a, b, c, d, e, f, g = os.path.basename(filename)[:-8].split('_')
        fileinfo = Fileinfo(a+b+c, d, int(e), int(f), int(g))

        injection = _calculate_current(fileinfo, IV, IF, IVseries)

        return cls(filename, fileinfo, injection, time, data, params)

    @property
    def time(self):
        return self.wave.x[-1]

    @property
    @utilities.once
    def baseline(self):
        return _find_baseline(self.wave,
                              before=self.params.baseline_before,
                              after=self.params.baseline_after)

    @property
    @utilities.once
    def steady(self):
        return _find_steady_state(self.wave,
                                  after=self.params.steady_after,
                                  before=self.params.steady_before,
                                  cutoff_percentile=self.params.steady_cutoff)

    @property
    @utilities.once
    def response(self):
        return self.steady - self.baseline

    @property
    @utilities.once
    def falling_curve(self):
        return _find_falling_curve(self.wave,
                                   window=self.params.falling_curve_window,
                                   before=self.params.steady_before)

    @property
    @utilities.once
    def falling_curve_fit(self):
        return _fit_falling_curve(self.falling_curve, self.baseline, self.steady)

    @property
    @utilities.once
    def rectification(self):
        return _find_rectification(self.falling_curve,
                                   self.steady,
                                   window_len=self.params.rectification_window)

    @property
    @utilities.once
    def _spike_i_and_threshold(self):
        return _find_spikes(self.wave)

    @property
    @utilities.once
    def spikes(self):
        spike_i, thresholds = self._spike_i_and_threshold
        return self.wave[spike_i]

    @property
    def spike_thresholds(self):
        spike_i, thresholds = self._spike_i_and_threshold
        return thresholds

    @property
    def spike_count(self):
        return len(self.spikes)

    @property
    @utilities.once
    def charging_curve_halfheight(self):
        "The height in the middle between depolarization and first spike"
        if self.spike_count < 1:
            return np.nan
        else:
            what = self.wave[(self.wave.x > self.params.steady_after)
                             & (self.wave.x < self.spikes[0].x)]
            return np.median(what.y)

    @property
    def depolarization_interval(self):
        return self.params.steady_before - self.params.steady_after

    @property
    @utilities.once
    def mean_isi(self):
        if self.spike_count > 2:
            return array_mean(np.diff(self.spikes.x))
        elif self.spike_count == 2:
            d = self.spikes.x[1]-self.spikes.x[0]
            return vartype(d, 0.001)
        else:
            return vartype(self.depolarization_interval, 0.001)

    @property
    @utilities.once
    def isi_spread(self):
        spikes = self.spikes
        if len(spikes) > 2:
            diff = np.diff(self.spikes.x)
            return diff.ptp()
        else:
            return np.nan

    @property
    @utilities.once
    def spike_latency(self):
        "Latency until the first spike or np.nan if no spikes"
        if len(self.spikes) > 0:
            return self.spikes[0].x-self.params.baseline_before
        else:
            return np.nan

    @property
    @utilities.once
    def mean_spike_height(self):
        return array_mean(self.spikes.y)

    @property
    @utilities.once
    def spike_height(self):
        spikes, thresholds = self._spike_i_and_threshold
        height = (self.spikes.y - thresholds)
        return height

    @property
    @utilities.once
    def spike_width(self):
        spikes, thresholds = self._spike_i_and_threshold
        ans = np.empty_like(spikes, dtype=float)
        x = self.wave.x
        y = self.wave.y
        halfheight = (self.spikes.y - thresholds) / 2 + thresholds
        for i, k in enumerate(spikes):
            beg = end = k
            while beg > 1 and y[beg - 1] > halfheight[i]:
                beg -= 1
            while end + 2 < y.size and y[end + 1] > halfheight[i]:
                end += 1
            ans[i] = (x[end] + x[end+1] - x[beg] - x[beg-1]) / 2
        return ans

    @property
    @utilities.once
    def spike_ahp(self):
        #do we need to use the spike threshold here also?
        spikes = self.spikes
        widths = self.spike_width
        x = self.wave.x
        y = self.wave.y
        ans = np.empty_like(spikes, dtype=float)
        for i in range(len(spikes)):
            #end should be the next spike (or when current inject stops)
            #begin should be close to spike peak location
            #1. find minimum 
            #2. AHP=MIN-THRESH
            end = spikes[i + 1].x if i < len(spikes)-1 else self.params.steady_before
            beg = min(spikes[i].x + widths[i], end)
            ans[i] = y[(x >= beg) & (x <= end)].min()
        return self.spike_thresholds - ans

class Attributable(object):
    _MEAN_ATTRIBUTES = {'mean_baseline', 'mean_spike_height', 'mean_spike_width', 'mean_spike_ahp'}
    _VAR_ARRAY_ATTRIBUTES = {'baseline', 'steady', 'response', 'rectification',
                             'mean_isi', 'spike_height'}
    _ARRAY_ATTRIBUTES = {'filename', 'injection', 'spike_width', 'spike_height',
                         'spike_latency', 'spike_count', 'spike_ahp',
                         'falling_curve_fit|params|amp', 'falling_curve_fit|params|tau',
                         'falling_curve_fit|function',
                         'charging_curve_halfheight',
                         'isi_spread'}

    def __getattribute__(self, attr):
        if attr == 'mean_spike_height':
            spikes = np.hstack([w.spikes for w in self])
            return array_mean(spikes['y'])
        elif attr == 'mean_spike_width':
            widths = np.hstack([w.spike_width for w in self])
            return array_mean(widths)
        elif attr == 'mean_spike_ahp':
            ahps = np.hstack([w.spike_ahp for w in self])
            return array_mean(ahps)
        elif attr in Attributable._MEAN_ATTRIBUTES:
            return vartype.average(getattr(self, attr[5:]))
        elif attr in Attributable._VAR_ARRAY_ATTRIBUTES:
            return vartype.array([getattr(wave, attr) for wave in self.waves])
        elif attr in Attributable._ARRAY_ATTRIBUTES:
            op = operator.attrgetter(attr.replace('|', '.'))
            ans = [op(wave) for wave in self.waves]
            return np.array(ans)
        else:
            return super(Attributable, self).__getattribute__(attr)

    def __getitem__(self, index):
        if isinstance(index, (slice, np.ndarray)):
            c = copy.copy(self)
            c.waves = self.waves[index]
            return c
        else:
            return self.waves[index]

class Measurement(Attributable):
#These IV and IF initial values may need to be adjusted, depending on the experiment.
#time is the duration of your IV or IF curve, need to adjust based on your data
    def __init__(self, dirname,
                 IV=(-500e-12, 50e-12),
                 IF=(200e-12, 20e-12),
                 time=.9,
                 bad_sweep=(),IVseries=4):
        self.name = os.path.basename(dirname)
        self.params = Params()

        ls = os.listdir(dirname)

        #calls IVCurve.load once for each file in dirname
        #thus, you need to have directories with only data in them
        #each experiment should be in separate directory, but IV and IF curves can be together
        waves = [IVCurve.load(dirname, f, IV, IF, self.params, IVseries,time=time)
                 for f in ls]
        self.waves = np.array([wave for wave in waves
                               if [wave.fileinfo.sweep,wave.fileinfo.series] not in bad_sweep])
        self.waves = self.waves[self.injection.argsort()]
