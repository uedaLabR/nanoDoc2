import os
import pandas as pd
from multiprocessing import Pool
from functools import partial
from nanoDoc2.utils import FileIO
from numba import jit, njit

binsize = 1000
def getBinkey(read):

    chrom = read.chrom
    start = read.r_st
    bin = start // binsize
    strand = read.strand
    #
    binkey = str(chrom)+"_"+str(strand)+"_"+str(bin)
    return binkey

def _writeToFile(item,pathout,filename):

    binkey,datalist = item
    os.makedirs(pathout + "/" + binkey, exist_ok=True)
    file_out = pathout + "/" + binkey + "/" + filename +"_pre.pq"
    df = pd.DataFrame(datalist,
                      columns=['chr', 'strand', 'start', 'end','cigar','fastq','offset','traceintervals','trace','signal'])

    #pd.to_pickle(df, file_out)
    FileIO.writeToPq(df, file_out)

import numpy as np
def boundaryToIntervals(traceboundary):

    arr = np.array(traceboundary)
    diff = np.diff(arr)
    # print(arr)
    # print("------")
    # print(diff)
    # print("======")
    diff = np.clip(diff,0,65535)
    return diff


def byteToHalf(b):
    return int((b / 256.0) * 16.0)

def to16bit(a_trace):


    _a = max(a_trace[0], a_trace[4])
    _c = max(a_trace[1], a_trace[5])
    _g = max(a_trace[2], a_trace[6])
    _u = max(a_trace[3], a_trace[7])
    _a = byteToHalf(_a)
    _c = byteToHalf(_c)
    _g = byteToHalf(_g)
    _u = byteToHalf(_u)
    #print(_a,_c,_g,_u)

    a_bi = _a << 12
    c_bi = _c << 8
    g_bi = _g << 4
    trace_bi = a_bi | c_bi | g_bi | _u
    #print(bin(trace_bi))
    return trace_bi

def convertTo16bit(trace):

    ret = []
    for t in trace:
        r = to16bit(t)
        ret.append(r)
    return ret


def toTuple(read):

    #
    offset = read.traceboundary[0]
    if offset < 0:
        offset = 0

    traceintervals = boundaryToIntervals(read.traceboundary)
    strand = read.strand
    if strand == -1:
        strand=0
    trace = convertTo16bit(read.trace)
    return read.chrom, strand, read.r_st, read.r_en,read.cigar_str,read.fastq,offset,traceintervals,trace,read.signal


def writeToFile(pathout,ncore,reads,filename):

    datadict = {}
    for read in reads:
        binkey = getBinkey(read)
        if binkey not in datadict:
            datadict[binkey] = []
        datadict[binkey].append(toTuple(read))

    filewriteF = partial(_writeToFile,pathout=pathout,filename=filename)
    with Pool(ncore) as p:
        p.map(filewriteF, datadict.items())



