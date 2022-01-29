import h5py
from statsmodels import robust
import numpy as np
import random
from scipy import interpolate
import pyarrow as pa
import pyarrow.parquet
import pandas as pd
from statistics import mean
import numba

def getMeans(signal,traceboundary):

    means = []
    tbb4 = 0
    unit = 10
    for tb in traceboundary:

        if (tb == 0) or (tbb4 == tb):
            continue
        if tb*unit >= len(signal):
            break
        subsignal = signal[tbb4*unit:tb*unit]
        if len(subsignal) >=1:
            means.append(mean(subsignal))
        tbb4 = tb

    return means

def theoryMean(fmercurrent,lgenome,strand):

    a = {}
    with open(fmercurrent) as f:
        cnt = 0
        for line in f:
            if cnt > 0:
                data = line.split()
                a[data[0]] = float(data[1])
            cnt = cnt+1

    means = []
    if strand == "-":

        rg = lgenome
        for n in range(0, len(rg) - 5):
            fmer = rg[n:n + 5]
            cv = a[fmer]
            means.append(cv)

    else:

        #plus strand
        rg = lgenome[::-1]
        for n in range(0,len(rg)-5):

           fmer = rg[n:n+5]
           cv = a[fmer]
           means.append(cv)

        means.reverse()

    return means

def predictShift(a,b):

    maxn = 0
    max = 0
    ar = None
    br = None
    for n in range(-3,3):

        bmod = moda(b,n)
        if len(a) > len(bmod):
            a = a[0:len(bmod)]
        if len(bmod) > len(a):
            bmod = bmod[0:len(a)]

        v = np.dot(a,bmod)
        if v > max:
            max = v
            maxn = n
            ar = a
            br = bmod

    return maxn,ar,br

def calcNormalizeScaleLMS(signalmeans,theorymean):

    try:
        y1 = 0
        y2 = 0
        y2_pow2 = 0
        y1y2 = 0
        n_len = len(signalmeans)
        for n in range(len(signalmeans)):

            y1 = y1 + theorymean[n]
            y2 = y2 + signalmeans[n]
            y1y2 = y1y2 + y1*y2
            y2_pow2 = y2_pow2 + y2*y2

        a = [[y2_pow2, y2], [y2, n_len]]
        ainv = np.linalg.inv(a)
        b = np.array([[y1y2], [y1]])
        return np.dot(ainv, b)
    except:
        return None

def moda(a,shift):

    if shift ==0:
        return a
    elif shift > 0:
        return a[shift:]
    else:
        a = list(a)
        for n in range(abs(shift)):
            a.insert(0, 0.5)
        return np.array(a)


def downsample(array, npts):
    interpolated = interp1d(np.arange(len(array)), array, axis = 0, fill_value = 'extrapolate')
    downsampled = interpolated(np.linspace(0, len(array), npts))
    return downsampled

def getFunction(scaleshifts,traceboundary,window,step):

    cnt = 0
    unit = 10

    x = []
    y1 = []
    y2 = []
    for v in  scaleshifts:
        if v is None:
            cnt+=step
            continue
        a,b = v

        if len(traceboundary) <= cnt+window:
            break

        x0 = (traceboundary[cnt]*unit + traceboundary[cnt+window]*unit)/2
        cnt+=step
        x.append(x0)
        y1.append(a)
        y2.append(b)

    f1 = interpolate.interp1d(x, y1, kind="quadratic",fill_value="extrapolate")
    f2 = interpolate.interp1d(x, y2, kind="quadratic",fill_value="extrapolate")
    return f1,f2


def normalizeSignal(read,traceboundary,fmercurrent):

    try:
        ret = _normalizeSignal(read,traceboundary,fmercurrent)
        #print("success")
        return ret

    except:
        print("error in normalization, use old method to normalize")
    return  normalizeSignal_old(read,traceboundary,fmercurrent)

from scipy.interpolate import interp1d
def _normalizeSignal(read,traceboundary,fmercurrent):

    low_limit = 40
    high_limit = 160

    lgenome = read.refgenome
    signal = read.signal
    strand = read.strand

    signalmeans = getMeans(signal,traceboundary)
    theorymean = theoryMean(fmercurrent, lgenome ,strand)
    shift, signalmeans, theorymean = predictShift(signalmeans, theorymean)
    window = 40
    step = 10
    start = 0
    end = window
    scaleshifts = []
    a_for_max_min = []
    b_for_max_min = []
    while (end + step) < len(signalmeans):

        scaleshift = calcNormalizeScaleLMS(signalmeans[start:end], theorymean[start:end])
        if scaleshift is None:
            scaleshifts.append(None)
        # by each 25 base 50 nt interval calculate a,b
        start = start + step
        end = end + step
        a = scaleshift[0][0]
        b = scaleshift[1][0]
        scaleshifts.append((a,b))
        a_for_max_min.append(a)
        b_for_max_min.append(b)


    functionA,functionB = getFunction(scaleshifts,traceboundary,window,step)
    num = np.arange(len(signal))
    a_ary = np.frompyfunc(functionA, 1, 1)(num) #mean
    b_ary = np.frompyfunc(functionB, 1, 1)(num) #sd

    # in order to avoid at adjustment with overfitting on both end
    a_ary = np.clip(a_ary, max(a_for_max_min), min(a_for_max_min))
    b_ary = np.clip(b_ary, max(b_for_max_min), min(b_for_max_min))

    signal = signal * a_ary + b_ary

    signal = np.clip(signal, low_limit, high_limit)
    signal = ((signal-low_limit) / (high_limit-low_limit)) * 255
    signal = np.around(signal.astype(np.float), 0)

    downsamplesize = len(signal) // 2 #half the size
    signal = downsample(signal, downsamplesize)

    signal = np.clip(signal, 0, 255)
    signal = signal.astype(np.uint8)
    return signal

def normalizeSignal_old(read,traceboundary,fmercurrent):

    low_limit = 50
    high_limit = 160

    lgenome = read.refgenome
    signal = read.signal
    strand = read.strand

    signalmeans = getMeans(signal,traceboundary)
    theorymean = theoryMean(fmercurrent, lgenome,strand)
    shift, signalmeans, theorymean = predictShift(signalmeans, theorymean)
    scaleshift = calcNormalizeScaleLMS(signalmeans, theorymean)
    #do it by 50 bp bin

    a = scaleshift[0][0]
    b = scaleshift[1][0]
    signal = signal * a + b
    signal = np.clip(signal, low_limit, high_limit)
    signal = signal - low_limit
    signal = (signal / (high_limit-low_limit)) * 255
    downsamplesize = len(signal) // 2 #half the size
    signal = downsample(signal, downsamplesize)
    signal = signal.astype(np.uint8)
    return signal
