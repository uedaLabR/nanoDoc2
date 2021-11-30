import glob
from sklearn.model_selection import train_test_split
import pyarrow.parquet as pq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import tensorflow as tf  # add
import numpy as np
import itertools
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import LocalOutlierFactor
from sklearn import metrics
from sklearn.preprocessing import MinMaxScaler
from Bio import SeqIO
from numpy import mean, absolute
import math
import gc
from scipy.optimize import curve_fit
import os
from numba import jit,u1,i8,f8

DATA_LENGTH_UNIT = 60
DATA_LENGTH = DATA_LENGTH_UNIT * 5 + 20


@jit(nopython=True)
def noise(sigma):
    return random.gauss(0, sigma)


from scipy import ndimage as ndi
import statistics


def binned(trimsignal, trimlength, mode=1):

    if len(trimsignal) == trimlength:
        return trimsignal  # not very likely to happen

    if len(trimsignal) > trimlength:
        # trim from first
        lefthalf = (len(trimsignal) - trimlength) // 2 - 1
        return trimsignal[lefthalf:trimlength]
    else:
        #
        ret = np.zeros(trimlength)
        diff = np.array([1, 0, -1])
        trimsignal = trimsignal.astype(np.float32)
        med = statistics.median(trimsignal)
        diffsig = ndi.convolve(trimsignal, diff)
        sigma = np.std(diffsig) / 10

        siglen = len(trimsignal)
        left = trimlength - siglen
        lefthalf = left // 2
        #
        for n in range(trimlength):
            if n < lefthalf or n >= (trimlength - lefthalf):
                if mode == 1:
                    ret[n] = med + noise(sigma)
                else:
                    ret[n] = 0  # zero pad
            else:
                idxa = n - lefthalf - 1
                if idxa < 0:
                    idxa = 0
                ret[n] = trimsignal[idxa]
        return ret

def binnedtrace(trace, sgstarttrace, sgendtrace, binsize):

    ret = []
    for t in trace:

        print(t,len(t),sgstarttrace,sgendtrace)
        if len(t) > 0:
            t = t[sgstarttrace:sgendtrace]
            print(t)
            t = binned(t, binsize, mode=0)
        ret.append(t)

    return ret


class PqReader:

    def __init__(self, indexfile, minreadlen):

        tb = []
        self.indexfile = indexfile
        self.minreadlen = minreadlen
        self.pq = None
        indexfile = indexfile + "/index.txt"
        # print(indexfile)
        f = open(indexfile)
        idxno = 0
        rowcnt = 0
        chrb4 = ""
        for s_line in f:
            s_line = s_line.split(",")

            if chrb4 != s_line[1]:
                idxno = idxno + 1
                rowcnt = 0

            chrb4 = s_line[1]

            rowcnt = rowcnt + 1
            if (rowcnt % 4000 == 0):
                idxno = idxno + 1

            tb.append((idxno, s_line[1], int(s_line[2]), int(s_line[3])))

        f.close()
        self.df = pd.DataFrame(tb, columns=('idx', 'chr', 'start', 'end'))
        # print(self.df)
        self.currentPqs = set()

    def loadpq(self, set0, chr, pos, strand):

        ss = sorted(set0)
        cnt = 0
        if len(self.currentPqs) > 0 and self.currentPqs <= set0:
            ss = set0 - self.currentPqs
            cnt = len(self.currentPqs)
            totaldf = self.pq

        totaldf = None
        for s in ss:

            s1 = self.indexfile + "/algined" + str(s) + ".pq"
            if not os.path.exists(s1):
                s1 = self.indexfile + "/" + str(s) + ".pq"

            print(s1)
            # columns = ['mapped_chrom', 'mapped_strand', 'mapped_start', 'mapped_end', 'clipped_bases_start',
            #            'clipped_bases_end', 'fastq', 'chrom', 'st', 'ed', 'cigar', 'md', 'move', 'trace', 'signal',
            #            'sgmean', 'sgst', 'sglen', 'bases'
            #            ])
            table = pq.read_table(s1,
                                  columns=['mapped_chrom', 'mapped_strand', 'mapped_start',
                                           'mapped_end', 'signal', 'trace','sgmean','sgst','sglen'])
            df = table.to_pandas()
            if strand == "+":
                df = df.query('(mapped_end - 40) >= ' + str(
                    pos) + ' & mapped_chrom=="' + chr + '" & mapped_strand =="' + strand + '" & (mapped_end - mapped_start) > ' + str(
                    self.minreadlen))
            else:
                df = df.query('(mapped_start +40) <= ' + str(
                    pos) + ' & mapped_chrom=="' + chr + '" & mapped_strand =="' + strand + '" & (mapped_end - mapped_start) > ' + str(
                    self.minreadlen))
            gc.collect()
            if cnt == 0:
                totaldf = df
            else:
                totaldf = pd.concat([totaldf, df], axis=0)

            cnt = cnt + 1

        return totaldf

    def freeUnsued(self, chr, pos, strand):

        if strand == "+":
            df = self.df.query('(mapped_end - 40) >= ' + str(
                pos) + ' & mapped_chrom=="' + chr + '" & mapped_strand =="' + strand + '" & (mapped_end - mapped_start) > ' + str(
                self.minreadlen))
        else:
            df = self.df.query('(mapped_start +40) <= ' + str(
                pos) + ' & mapped_chrom=="' + chr + '" & mapped_strand =="' + strand + '" & (mapped_end - mapped_start) > ' + str(
                self.minreadlen))

    def checkUpdate(self, chr, pos, strand):

        df = self.df.query('start <= ' + str(pos) + ' & end >= ' + str(pos + 5) + ' & chr == "' + chr + '"')
        df = df['idx'].drop_duplicates()
        s = set(df)
        if (s != self.currentPqs or self.pq is None):
            self.pq = self.loadpq(s, chr, pos, strand)
            self.currentPqs = s

        # if pos%10 ==0:
        #    self.freeUnsued(chr,pos,strand)

    #             print(pos,s)


    def getData(self, chr, strand, pos, maxtake):

        print('maxtake',maxtake,chr,pos)
        unitwidth = DATA_LENGTH_UNIT
        self.checkUpdate(chr, pos, strand)
        if self.pq is None:
            return None

        pqlen = len(self.pq)
        _pq = self.pq
        # print("pqlen=",len(_pq))
        data = []
        tracedata = []
        takecnt = 0
        for index, row in _pq.iterrows():

            mapped_chrom = row['mapped_chrom']
            mapped_start = row['mapped_start']
            mapped_end = row['mapped_end']
            #'signal', 'trace', 'sgmean', 'sgst', 'sglen'
            sgst = row['sgst']

            if strand == "-":

                relpos = mapped_end - pos
                if relpos < 0:
                    continue
                if pos - 7 < mapped_start:
                    continue
                if relpos+7 >= len(sgst):
                    continue

                sgstart = sgst[relpos]
                sgend = sgst[relpos+7]
                signal = row['signal'][sgstart:sgend]
                if signal is not None and len(signal)> 50:
                    signal = binned(signal, 420)
                if signal is not None and len(signal) == 420:
                    data.append(signal)
                    takecnt = takecnt + 1

                # trace = row['trace']
                # #
                # trace = np.reshape(trace, (8, -1))
                # #print(trace.shape)
                # sgstarttrace = (sgstart - sgst[0]) // 10
                # sgendtrace = (sgend - sgst[0]) // 10
                #
                # trace = binnedtrace(trace, sgstarttrace, sgendtrace, 42)
                # tracedata.append(trace)

            else:

                relpos = pos - mapped_start
                if relpos < 0:
                    continue
                if pos + 7 > mapped_end:
                    continue
                if relpos+7 >= len(sgst):
                    continue

                sgstart = sgst[relpos]
                sgend = sgst[relpos+7]
                signal = row['signal'][sgstart:sgend]
                if signal is not None and len(signal)> 50:
                    signal = binned(signal, 420)
                if signal is not None and len(signal) == 420:
                    data.append((takecnt,signal))
                    takecnt = takecnt + 1

                # trace = row['trace']
                # #
                # trace = np.reshape(trace, (8, -1))
                # #print(trace.shape)
                #
                # sgstarttrace = (sgstart - sgst[0]) // 10
                # sgendtrace = (sgend - sgst[0]) // 10
                # trace = binnedtrace(trace,sgstarttrace,sgendtrace,42)
                # tracedata.append(trace)

            if takecnt == maxtake:
                break

        print(chr,pos,len(data),maxtake)
        return data,len(data)


