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
import pysam

from Bio import SeqIO
from numpy import mean, absolute
import math
import gc
from scipy.optimize import curve_fit
import os
from numba import jit,u1,i8,f8

DATA_LENGTH_UNIT = 30
DATA_LENGTH = DATA_LENGTH_UNIT * 5 + 10
TraceToSignalRatio = 5
binsize = 1000
takemargin = 10

from scipy import ndimage as ndi
import statistics


def binned(trimsignal, trimlength, mode=1):

    if len(trimsignal) == trimlength:
        return trimsignal  # not very likely to happen

    if len(trimsignal) > trimlength:
        # trim from first
        lefthalf = ((len(trimsignal) - trimlength)) // 2 - 1
        #print("trim",len(trimsignal),trimlength,lefthalf,len(trimsignal[lefthalf:lefthalf+trimlength]))
        return trimsignal[lefthalf:lefthalf+trimlength]

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

def intervalToAbsolute(intervals):

    ret = []
    cnt = 0
    sum = 0
    for n in intervals:

        if cnt==0:
            ret.append(0)
            sum = sum + n
            cnt += 1
            continue
        #else
        ret.append(sum)
        sum = sum + n

    ret.append(sum)
    return np.array(ret)

import pyarrow.parquet as pq

class PqReader:

    def __init__(self, path,maxreads = 2000):

        self.path = path
        self.batch = None
        self.maxreads = int(maxreads * 1.2)# take little more sample since some reads disqualify
        self.maxreads_org = maxreads
        self.bufferData = None
        self.loadchr = None
        sortedfile = sorted(glob.glob(path+"/*.pq"))
        #make index from parquet meta data
        indexlist = []
        fileidx = 0
        for file in sortedfile:

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
            chrInfo = parquet_file.metadata.row_group(0).column(2).statistics.min
            strandInfo = parquet_file.metadata.row_group(0).column(3).statistics.min
            startInfo = parquet_file.metadata.row_group(0).column(4).statistics.min
            endInfo = parquet_file.metadata.row_group(0).column(5).statistics.max
            indexlist.append((fileidx,chrInfo,strandInfo,startInfo,endInfo))
            fileidx +=1

        self.indexdf = pd.DataFrame(indexlist,
                          columns=['fileidx','chr','strand','start','end'])

    def getDepth(self,chr, pos, strand):

        #extract parquet file contain reads in this region
        query = 'start <= ' + str(pos) + ' & end >= ' + str(pos) + ' & chr == "' + chr + '" & strand == ' + str(
            strand) + ''
        pqfiles = self.indexdf.query(query)
        sortedfile = sorted(glob.glob(self.path + "/*.pq"))
        indexdata = None
        for index, row in pqfiles.iterrows():

            fileidx = row['fileidx']
            filepath = sortedfile[fileidx]

            if indexdata is None:
                indexdata = pq.read_table(filepath, columns=['start', 'end']).to_pandas()


            else:
                dataadd = pq.read_table(filepath, columns=['start', 'end']).to_pandas()
                indexdata = pd.concat([indexdata, dataadd])

        depth = indexdata.query('start <=' + str(pos-10) + ' & end >=' + str(pos+10))['start'].count()
        return depth



    def randomsample(self, pos, data, ntake, indexes):

        if indexes is None:

            #print('init read')
            datainpos = data.query('start <=' + str(pos-takemargin) + ' & end >=' + str(pos+takemargin))
            if datainpos is None or len(datainpos) ==0:
                return None
            datainpos = datainpos.sample(n=ntake)
            return datainpos.loc[:, ['fileidx', 'read_no']]

        else:

            datainpos = data.query('start <=' + str(pos-takemargin) + ' & end >=' + str(pos+takemargin))
            df_alreadyhave = indexes['read_no']
            datainpos = datainpos[~datainpos.read_no.isin(df_alreadyhave)]
            cnt = ntake - len(df_alreadyhave)
            if (cnt > 0) and (cnt <= len(datainpos)):
                datainpos = datainpos.sample(n=cnt)
                #print(cnt,len(datainpos))
                return datainpos.loc[:, ['fileidx', 'read_no']]
            else:
                #print('return none')
                return None



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

        self.loadchr = chr
        #extract parquet file contain reads in this region
        query = 'start <= ' + str(pos) + ' & end >= ' + str(pos) + ' & chr == "' + chr + '" & strand == ' + str(
            strand) + ''
        pqfiles = self.indexdf.query(query)
        #
        sortedfile = sorted(glob.glob(self.path + "/*.pq"))
        #
        indexdata = None
        for index, row in pqfiles.iterrows():

            fileidx = row['fileidx']
            filepath = sortedfile[fileidx]

            if indexdata is None:
                indexdata = pq.read_table(filepath, columns=['read_no', 'chr', 'strand', 'start', 'end']).to_pandas()
                indexdata['fileidx'] = fileidx

            else:
                dataadd = pq.read_table(filepath, columns=['read_no', 'chr', 'strand', 'start', 'end']).to_pandas()
                dataadd['fileidx'] = fileidx
                indexdata = pd.concat([indexdata, dataadd])

        # get reads ids for bufferted position (start,end)
        start = pos
        end = ((pos // binsize)+1) * binsize
        if not strand:
            start = (pos // binsize) * binsize
            end = pos

        readsIndex = None
        ntake = self.maxreads
        # get readid to reads for bin interval batch
        for n in range(start,end):

            addIndex = self.randomsample(pos, indexdata, ntake, readsIndex)
            if readsIndex is None:
                readsIndex = addIndex
            elif addIndex is not None:
                readsIndex = pd.concat([readsIndex, addIndex])

        #print("start reading file")
        # read data with row signal
        fileindexes = readsIndex['fileidx'].unique()
        dataWithRow = None
        for findex in fileindexes:

            filepath = sortedfile[findex]
            readids = readsIndex['read_no']
            filterlist = []
            filterTp = ('read_no', 'in', readids)
            filterlist.append(filterTp)
            columns = ['read_no', 'chr', 'strand', 'start', 'end', 'cigar', 'offset', 'traceintervals', 'signal']
            if dataWithRow is None:
                dataWithRow = pq.read_table(filepath, filters=filterlist, columns=columns).to_pandas()

            else:
                dataadd = pq.read_table(filepath, filters=filterlist, columns=columns).to_pandas()
                dataWithRow = pd.concat([dataWithRow, dataadd])

        dataWithRow['traceintervals'] = dataWithRow['traceintervals'] .apply(intervalToAbsolute)
        self.bufferData = dataWithRow



    def getRowData(self, chr, strand, pos,takecnt=-1):

        if ((self.bufferData is None) or (chr != self.loadchr) or ((pos % binsize) == 0)):
            #print("loading row files",self.bufferData is None, chr != self.loadchr , (pos % binsize) == 0)
            self.load(chr, pos, strand)

        sampled,sampledlen = self.getFormattedData(strand, pos,takecnt)
        if sampledlen > self.maxreads_org:
            sampled = sampled[0:self.maxreads_org*DATA_LENGTH]
            sampledlen = self.maxreads_org
        return sampled,sampledlen

    import pysam
    def correctCigar(self,targetPos,cigar):

        a = pysam.AlignedSegment()
        a.cigarstring = cigar
        refpos = 0
        relpos = 0
        for cigaroprator, cigarlen in a.cigar:

            if cigaroprator == 0:  # match

                if refpos + cigarlen > targetPos:
                    return relpos + (targetPos - refpos)

                relpos = relpos + cigarlen
                refpos = refpos + cigarlen

            elif cigaroprator == 2:  # Del
                refpos = refpos + cigarlen
            elif cigaroprator == 1 or cigaroprator == 3:  # Ins or N
                relpos = relpos + cigarlen

        return 0

    def getRelativePos(self,strand,start,end,cigar,pos,traceintervalLen):

        #tp = (strand,start,end,cigar,pos,traceintervalLen)
        rel = pos - start
        rel = self.correctCigar(rel,cigar)
        start = rel
        if start < 2:
            # do not use lower end
            return None

        end = rel + 5
        if end >= traceintervalLen:
            end = traceintervalLen-1
        if start == end:
            return None

        # print(tp)
        # print("start end",start,end)
        return start,end


    def calcStartEnd(self,strand,start,end,cigar,pos,offset,traceintervals):

        rp = self.getRelativePos(strand,start,end,cigar,pos,len(traceintervals))
        if rp is None:
            return None

        relativeStart, relativeEnd = rp
        # print("relativeStart, relativeEnd",relativeStart, relativeEnd)
        # print(traceintervals)
        if relativeEnd >= len(traceintervals) or relativeStart >= len(traceintervals):
            return None
        sig_start = offset + traceintervals[relativeStart] * TraceToSignalRatio
        sig_end = offset + traceintervals[relativeEnd] * TraceToSignalRatio
        return sig_start,sig_end


    def getOneRow(self,row,strand, pos):

        chr = row['chr']
        start = row['start']
        end = row['end']
        cigar = row['cigar']
        offset = row['offset']
        traceintervals = row['traceintervals']
        signal = row['signal']
        #
        sted = self.calcStartEnd(strand,start,end,cigar,pos,offset,traceintervals)
        if sted is None:
            return None
        sig_start, sig_end = sted
        #print(sig_start, sig_end)
        signal = signal[sig_start:sig_end]
        #print(signal)
        #
        binnedSignal = binned(signal, DATA_LENGTH)
        #print('bin signal len',len(binnedSignal))
        return binnedSignal

    def getFormattedData(self, strand, pos,_takecnt):

        data = []
        takecnt = 0
        for index, row in self.bufferData.iterrows():

            rowdata = self.getOneRow(row, strand, pos)
            if rowdata  is not None:
                data.extend(rowdata)
                takecnt = takecnt + 1
                #print("takecnt",takecnt,_takecnt)
                if _takecnt > 0 and takecnt == _takecnt:
                    break

        data = np.array(data)
        data = data / 256
        return  data,takecnt
