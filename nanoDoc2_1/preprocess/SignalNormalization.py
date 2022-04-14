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


import pysam
def getrelpos(cigar,pos):
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    refpos = 0
    relpos = 0
    for cigaroprator, cigarlen in a.cigar:

        if cigaroprator == 0:  # match

            if refpos + cigarlen > pos:
                return relpos + (pos - refpos)

            relpos = relpos + cigarlen
            refpos = refpos + cigarlen

        elif cigaroprator == 2:  # Del
            refpos = refpos + cigarlen
        elif cigaroprator == 1 or cigaroprator == 3:  # Ins or N
            relpos = relpos + cigarlen

    return 0

from statistics import mean
def getMeans(signal,traceboundary,cigar,seqlen):

    means = []
    unit = 10
    for n in range(0, seqlen - 5):

        relpos = getrelpos(cigar,n)
        if relpos+1 < len(traceboundary):
            start = traceboundary[relpos] * unit
            end = traceboundary[relpos+1] * unit
            if end < len(signal):
                subsignal = signal[start:end]
                if len(subsignal) > 0:
                    means.append(mean(subsignal))


    return means


def theoryMean(fmerDict,lgenome,strand):

    means = []
    if strand == "-":

        rg = lgenome
        for n in range(0, len(rg) - 5):
            fmer = rg[n:n + 5]
            if "N" in fmer:
                fmer = fmer.replace('N', 'A')
            cv = fmerDict[fmer]
            means.append(cv)

    else:

        #plus strand
        rg = lgenome[::-1]
        for n in range(0,len(rg)-5):

           fmer = rg[n:n+5]
           if "N" in fmer:
               fmer = fmer.replace('N', 'A')
           cv = fmerDict[fmer]
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

from numba import jit

def average(arr, n):
    end =  n * int(len(arr)/n)
    return np.mean(arr[:end].reshape(-1, n), 1).astype('float32')

import scipy.signal as scisignal
def downsample(array, npts):

    # interpolated = interp1d(np.arange(len(array)), array, axis = 0, fill_value = 'extrapolate')
    # downsampled = interpolated(np.linspace(0, len(array), npts))

    #downsampled = scisignal.resample(array, npts)
    downsampled = average(array, 2)
    return downsampled



def getFunction(scaleshifts,traceboundary,window,step):

    cnt = 0
    unit = 5 # apply after down sampled

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

    try:
        f1 = interpolate.interp1d(x, y1, kind="quadratic",fill_value="extrapolate")
        f2 = interpolate.interp1d(x, y2, kind="quadratic",fill_value="extrapolate")
        #print("pass1")

    except  Exception as e:

        f1 = interpolate.interp1d(x, y1, kind="slinear",fill_value="extrapolate")
        f2 = interpolate.interp1d(x, y2, kind="slinear",fill_value="extrapolate")
        #print("pass2")

    return f1,f2


def normalizeSignal(read,traceboundary,fmerDict):

    data = None
    try:
        data =  _normalizeSignal(read,traceboundary,fmerDict)
    except:

        print("normalize by window failed with " + read.read_id +" len="+str(len(read.sequence))+" , probably too short,sequence is normized as a whole")
        data = normalizeSignal_as_whole(read,traceboundary,fmerDict)

    #data = normalizeSignal_as_whole(read, traceboundary, fmerDict)
    # if data is not None:
    data = format(data)
    #data = format(data)
    return  data

def format(signal):

    low_limit = 40
    high_limit = 160


    signal = np.clip(signal, low_limit, high_limit)
    signal = ((signal-low_limit) / (high_limit-low_limit)) * 255
    signal = np.around(signal.astype(np.float), 0)

    signal = np.clip(signal, 0, 255)
    signal = signal.astype(np.uint8)

    return signal


def format(signal):

    low_limit = 40
    high_limit = 160

    #downsamplesize = len(signal) // 2 #half the size
    #signal = downsample(signal, downsamplesize)

    signal = np.clip(signal, low_limit, high_limit)
    signal = ((signal-low_limit) / (high_limit-low_limit)) * 255
    signal = np.around(signal, 0)
    signal = np.clip(signal, 0, 255)
    signal = signal.astype(np.uint8)

    return signal


from scipy.interpolate import interp1d
def _normalizeSignal(read,traceboundary,fmerDict):



    lgenome = read.refgenome
    signal = read.signal
    strand = read.strand
    cigar = read.cigar_str

    signalmeans = getMeans(signal,traceboundary,cigar,len(lgenome))

    theorymean = theoryMean(fmerDict, lgenome ,strand)
    shift, signalmeans, theorymean = predictShift(signalmeans, theorymean)
    if signalmeans is None:
        return signal

    window = 60
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
        # by each 5 base 40 nt interval calculate a,b
        start = start + step
        end = end + step
        a = scaleshift[0][0]
        b = scaleshift[1][0]
        scaleshifts.append((a,b))
        a_for_max_min.append(a)
        b_for_max_min.append(b)
    #


    functionA,functionB = getFunction(scaleshifts,traceboundary,window,step)
    num = np.arange(len(signal))
    a_ary = functionA(num) #mean
    b_ary = functionB(num) #sd

    # in order to avoid at adjustment with overfitting on both end
    a_ary = np.clip(a_ary, min(a_for_max_min), max(a_for_max_min))
    b_ary = np.clip(b_ary, min(b_for_max_min), max(b_for_max_min))
    signal = signal * a_ary + b_ary

    return signal

def normalizeSignal_as_whole(read,traceboundary,fmerDict):

    lgenome = read.refgenome
    signal = read.signal
    strand = read.strand
    cigar = read.cigar_str

    signalmeans = getMeans(signal,traceboundary,cigar,len(lgenome))
    theorymean = theoryMean(fmerDict, lgenome,strand)
    shift, signalmeans, theorymean = predictShift(signalmeans, theorymean)

    if signalmeans is None:
        return signal

    scaleshift = calcNormalizeScaleLMS(signalmeans, theorymean)

    if scaleshift is None:
        return signal

    a = scaleshift[0][0]
    b = scaleshift[1][0]
    signal = signal * a + b

    return signal
