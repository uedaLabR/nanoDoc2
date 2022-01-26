import glob

import pyarrow.parquet as pq
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import numpy as np
import itertools
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from numpy import mean, absolute
import math
import gc
from scipy.optimize import curve_fit
import os
from numba import jit,u1,i8,f8

DATA_LENGTH_UNIT = 60
DATA_LENGTH = DATA_LENGTH_UNIT * 5 + 20

binsize = 1000
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

        rand1 = np.random.rand(lefthalf) * sigma
        rand1 = rand1 + med
        leftlen = trimlength - siglen - lefthalf
        rand2 = np.random.rand(leftlen) * sigma
        rand2 = rand2 + med
        #
        ret = np.concatenate([rand1, trimsignal, rand2])

        return ret

@jit(nopython=True)
def binnedtrace(trace, sgstarttrace, sgendtrace, binsize):

    ret = []
    for t in trace:

        if len(t) > 0:
            t = t[sgstarttrace:sgendtrace]
            t = binned(t, binsize, mode=0)
        ret.append(t)

    return ret

import pyarrow.parquet as pq

class PqReader:

    def __init__(self, path,maxreads = 2000):

        self.path = path
        self.batch = None
        self.maxreads = maxreads
        self.currentData = None
        filelist = glob.glob(path+"/*.pq")
        #make index from parquet meta data
        self.df = None
        self.indexes = None
        indexlist = []
        for file in filelist:

            # ('read_id', string()),
            # ('chr', string()),
            # ('strand', bool_()),
            # ('start', uint32()),
            # ('end', uint32()),
            # ('cigar', string()),
            # ('fastq', string()),
            # ('offset', uint16()),
            # ('traceintervals', list_(uint16())),
            # ('trace', list_(uint16())),
            # ('signal', list_(uint8()))
            parquet_file = pq.ParquetFile(file)
            chrInfo = parquet_file.metadata.row_group(0).column(1).statistics.min
            strandInfo = parquet_file.metadata.row_group(0).column(2).statistics.min
            startInfo = parquet_file.metadata.row_group(0).column(3).statistics.min
            endInfo = parquet_file.metadata.row_group(0).column(4).statistics.max
            indexlist.append((file,chrInfo,strandInfo,startInfo,endInfo))

        self.df = pd.DataFrame(indexlist,
                          columns=['filepath', 'chr','strand','start','end'])



    # ('read_id', string()),
    # ('chr', string()),
    # ('strand', bool_()),
    # ('start', uint32()),
    # ('end', uint32()),
    # ('cigar', string()),
    # ('fastq', string()),
    # ('offset', uint16()),
    # ('traceintervals', list_(uint16())),
    # ('trace', list_(uint16())),
    # ('signal', list_(uint8()))
    def load(self,chr, pos, strand):

        self.alldata = None
        query = 'start <= ' + str(pos) + ' & end >= ' + str(pos) + ' & chr == "' + chr + '" & strand == ' + str(strand) + ''
        pqfiles = self.df.query(query)
        data = None
        for index, row in pqfiles.iterrows():

            filepath = row['filepath']
            print(filepath)
            if data is None:
                data = pq.read_table(filepath, columns=['read_id','chr', 'strand','start','end']).to_pandas()
                data["file"] = filepath
            else:
                dataadd = pq.read_table(filepath, columns=['read_id','chr', 'strand','start','end']).to_pandas()
                dataadd["file"] = filepath
                data = pd.concat([data,dataadd])

        self.indexdata = data
        self.indexes = self.randomsample(pos,data,self.maxreads)
        alldata = None
        for n, row in pqfiles.iterrows():

            filepath = row['filepath']
            readids = self.indexes.query('file == ' +"'"+filepath+"'")['read_id'].tolist()
            #print(readids)
            filterlist = []
            filterTp = ('read_id','in',readids)
            filterlist.append(filterTp)
            if alldata is None:
                alldata = pq.read_table(filepath, filters=filterlist, columns=['read_id', 'chr', 'start', 'end','cigar','offset','traceintervals','signal']).to_pandas()
            else:
                alldataadd = pq.read_table(filepath, filters=filterlist,
                                        columns=['read_id', 'chr', 'start', 'end', 'cigar', 'offset', 'traceintervals',
                                                 'signal']).to_pandas()
                alldata = pd.concat([alldata, alldataadd])

        self.currentData = alldata

    def loadAdditional(self,addindexes):


        pqfiles = addindexes['file'].unique()
        alldata = None
        for filepath in pqfiles:

            readids = addindexes.query('file == ' + "'" + filepath + "'")['read_id'].tolist()
            # print(readids)
            filterlist = []
            filterTp = ('read_id', 'in', readids)
            filterlist.append(filterTp)

            if alldata is None:
                alldata = pq.read_table(filepath, filters=filterlist,
                                        columns=['read_id', 'chr', 'start', 'end', 'cigar', 'offset', 'traceintervals',
                                                 'signal']).to_pandas()
            else:
                alldataadd = pq.read_table(filepath, filters=filterlist,
                                           columns=['read_id', 'chr', 'start', 'end', 'cigar', 'offset',
                                                    'traceintervals',
                                                    'signal']).to_pandas()
                alldata = pd.concat([alldata, alldataadd])

        if alldata is not None:
            self.currentData = pd.concat([self.currentData, alldata])



    def randomsample(self,pos,data,ntake):

        datainpos = data.query('start <='+ str(pos) +' & end >=' + str(pos))
        datainpos = datainpos.sample(n=ntake)
        return datainpos.loc[:,['file','read_id']]

    def additional_randomsample(self,pos,data,ntake,indexes):

        datainpos = data.query('start <='+ str(pos) +' & end >=' + str(pos))
        df_alreadyhave = indexes['read_id']
        datainpos = datainpos[~datainpos.read_id.isin(df_alreadyhave)]
        datainpos = datainpos.sample(n=ntake)
        return datainpos.loc[:,['file','read_id']]

    def checkUpdate(self,pos):

        #check depth
        count = self.currentData.query('start <='+ str(pos) +' & end >=' + str(pos))['read_id'].count()
        print(count)
        if count == self.maxreads:
            return
        countInIndex = self.indexdata.query('start <=' + str(pos) + ' & end >=' + str(pos))['read_id'].count()
        if count == countInIndex:
            return
        #update
        takecount = countInIndex - count
        #
        addindexes = self.additional_randomsample(pos, self.indexdata, takecount,self.indexes)
        self.loadAdditional(addindexes)

        self.indexes = pd.concat([self.indexes, addindexes])
        #




    def getRowData(self, chr, strand, pos):

        if self.currentData is None:
            self.load(chr, pos, strand)
        else:
            self.checkUpdate(pos)
        print(len(self.currentData))




    def getData(self, chr, strand, pos, maxtake):

        print('maxtake',maxtake,chr,pos)
        unitwidth = DATA_LENGTH_UNIT
        self.checkUpdate(chr, pos, strand)
        if self.pq is None:
            return None

        pqlen = len(self.pq)
        _pq = self.pq
        print("pqlen=",len(_pq))
        data = []
        tracedata = []
        takecnt = 0
        dftake = _pq.sample(frac=1)
        for index, row in dftake.iterrows():

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


