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

DATA_LENGTH_UNIT = 8
SEQ_TAKE_MARGIN = 3
binsize = 2000
takemargin = 10

from scipy import ndimage as ndi
import statistics
@jit
def eachScore(nuc,atrace):


    if nuc == 'A':
        nucidx = 0
    elif nuc == 'C':
        nucidx = 1
    elif nuc == 'G':
        nucidx = 2
    else:
        nucidx = 3

    rowsum =np.sum(atrace, axis=1)
    sum = np.sum(rowsum)
    matchsum = rowsum[nucidx]
    return matchsum/sum

def getScore(subtraces,seq):

    score = 0
    upto = min(len(seq),len(subtraces))
    for n in range(upto):
        nuc = seq[n]
        atrace = subtraces[n]
        ascore = eachScore(nuc,atrace)
        score +=ascore
    return score


def analyzeIdxShift(subtraces,rseq):

    maxidx = 0
    maxscore = 0
    for n in range(6):
        seq = rseq[n:n+10]
        score = getScore(subtraces,seq)
        idx = n-3
        if score > maxscore:
            maxidx = idx

    #print("maxidx",maxidx)
    return maxidx


def binned(trimtrace, traceItv,trimUnitLength,rseq):

    start = 0
    prev = 0
    subtraces = []
    for n in range(len(traceItv)):
        if n == 0:
            start = traceItv[0]
            prev = 0
            continue
        s = prev
        e =  traceItv[n] - start
        subtrace = trimtrace[::,s:e]
        prev = e
        subtraces.append(subtrace)

    #analizeIdx
    idx = analyzeIdxShift(subtraces,rseq)
    nret = None
    for m in range(idx,min(idx+5,len(subtraces))):

        if m>=0 and m < len(subtraces):
            subtrace = subtraces[m]
            bineach = binnedEach(subtrace, trimUnitLength)
            #print("bineach",bineach)
            if nret is None:
                nret = bineach
            else:
                nret = np.concatenate((nret, bineach), axis=1)

    #print("n",n,nret.shape[1])
    datalen = 0
    if nret is not None:
        datalen = nret.shape[1]
    leftlen = trimUnitLength * 5 - datalen
    if leftlen > 0:
        zeros = np.zeros((4, leftlen))
        if nret is not None:
            nret = np.concatenate((nret, zeros), axis=1)
        else:
            nret = zeros

    return nret


def binnedEach(trimtrace, trimlength):

    tracelen = trimtrace.shape[1]

    if tracelen == trimlength:
        return trimtrace  # not very likely to happen

    if tracelen > trimlength:
        # trim from first
        lefthalf = (tracelen - trimlength) // 2 - 1
        return trimtrace[::,lefthalf:lefthalf+trimlength]

    else:
        #

        left = trimlength - tracelen
        lefthalf = left // 2
        zeroleft = [0] * lefthalf
        leftlen = trimlength - tracelen - lefthalf
        zeroright = [0] * leftlen
        ret = []
        for row in trimtrace:
            data = np.concatenate([zeroleft, row, zeroright])
            ret.append(data)

        ret = np.array(ret)
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
import mappy as mp

class PqReader:

    def getFilePathList(self,path):

        sortedfile = sorted(glob.glob(path + "/*.pq"))
        # print(sortedfile)
        pluslist = []
        minuslist = []
        for f in sortedfile:
            if "_-1_" in f:
                minuslist.append(f)
            else:
                pluslist.append(f)
        pluslist.sort(key=lambda x: int(x.split('_')[-1].replace(".pq", "")))
        minuslist.sort(key=lambda x: int(x.split('_')[-1].replace(".pq", "")))
        sortedfile = pluslist
        sortedfile.extend(minuslist)
        return sortedfile



    def __init__(self, path,ref,maxreads = 1000,IndelStrict = False):

        self.IndelStrict = IndelStrict
        self.path = path
        self.batch = None
        print("ref",ref)
        self.a = mp.Aligner(ref)
        self.maxreads = int(maxreads * 1.2)# take little more sample since some reads disqualify
        self.maxreads_org = maxreads
        self.bufferData = None
        self.loadchr = None

        #make index from parquet meta data
        indexlist = []
        fileidx = 0
        sortedfile = self.getFilePathList(path)

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
        print(self.indexdf)

    def getDepth(self,chr, pos, strand):

        #extract parquet file contain reads in this region
        query = 'start <= ' + str(pos-takemargin) + ' & end >= ' + str(pos+takemargin) + ' & chr == "' + chr + '" & strand == ' + str(
            strand) + ''
        pqfiles = self.indexdf.query(query)
        sortedfile = self.getFilePathList(self.path)
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
            if ntake < len(datainpos):
                datainpos = datainpos.sample(n=ntake)
            return datainpos

        else:

            dataposprev = indexes.query('start <=' + str(pos-takemargin) + ' & end >=' + str(pos+takemargin))
            if len(dataposprev) == ntake:
                return None

            else:
                datainpos = data.query('start <=' + str(pos-takemargin) + ' & end >=' + str(pos+takemargin))
                df_alreadyhave = dataposprev['read_no']
                datainpos = datainpos[~datainpos.read_no.isin(df_alreadyhave)]
                cnt = ntake - len(df_alreadyhave)
                if (cnt > 0) and (cnt <= len(datainpos)):
                    datainpos = datainpos.sample(n=cnt)
                    return datainpos

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

        print("loding data")
        self.loadchr = chr
        #extract parquet file contain reads in this region
        query = 'start <= ' + str(pos-takemargin) + ' & end >= ' + str(pos+takemargin) + ' & chr == "' + chr + '" & strand == ' + str(
            strand) + ''
        pqfiles = self.indexdf.query(query)
        #
        sortedfile = self.getFilePathList(self.path)
        print(sortedfile)
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
        # need to buffer from equally from interval
        # middle = (start+end) // 2
        # l = list(range(start+1,middle))
        # l2 = list(range(middle,end-1))
        # ll = [middle,end,start]
        # ll.extend(l)
        # ll.extend(l2)

        for pos2 in range(start,end):

            addIndex = self.randomsample(pos2, indexdata, ntake, readsIndex)
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
            columns = ['read_no', 'chr', 'strand', 'start', 'end', 'cigar','genome','offset', 'traceintervals','trace']
            if dataWithRow is None:
                dataWithRow = pq.read_table(filepath, filters=filterlist, columns=columns).to_pandas()

            else:
                dataadd = pq.read_table(filepath, filters=filterlist, columns=columns).to_pandas()
                dataWithRow = pd.concat([dataWithRow, dataadd])

        dataWithRow['traceintervals'] = dataWithRow['traceintervals'] .apply(intervalToAbsolute)
        self.bufferData = dataWithRow


    def getRowData(self, chr, strand, pos,takecnt=-1):

        if strand == "+":
            strand = True
        if strand == "-":
            strand = False

        if pos - SEQ_TAKE_MARGIN < 0:
            return None
        ADDITONAL_MARGIN = 3
        margin = SEQ_TAKE_MARGIN+ADDITONAL_MARGIN
        rseq = self.a.seq(chr,start=(pos-margin),end=(pos+margin))

        if ((self.bufferData is None) or (chr != self.loadchr) or ((pos % binsize) == 0)):
            #print("loading row files",self.bufferData is None, chr != self.loadchr , (pos % binsize) == 0)
            print("load data init")
            self.bufferData = None
            self.load(chr, pos, strand)

        traces,sampledlen = self.getFormattedData(strand, pos,rseq,takecnt)
        if takecnt > self.maxreads_org:
            sampled = traces[0:self.maxreads_org]
            sampledlen = self.maxreads_org
            traces = sampled

        if sampledlen < takecnt and takecnt > 0:
            depth = self.getDepth(chr, pos, strand)
            if depth > sampledlen + 30:
                #load and sample again since not enough sampling for this region
                print("load data",depth,sampledlen)
                self.load(chr, pos, strand)
                traces,sampledlen = self.getFormattedData(strand, pos,rseq, takecnt)

        if takecnt == -1:
            depth = self.getDepth(chr, pos, strand)
            if sampledlen < depth:
                self.load(chr, pos, strand)
                traces,sampledlen = self.getFormattedData(strand, pos,rseq, takecnt)


        return traces,sampledlen


    def getRowSequence(self, chr, strand, pos,takecnt=-1):

        binsizehalf = binsize //2
        if ((self.bufferData is None) or (chr != self.loadchr) or ((pos % binsizehalf) == 0)):
            #print("loading row files",self.bufferData is None, chr != self.loadchr , (pos % binsize) == 0)
            self.load(chr, pos, strand)

        return self.bufferData.iterrows()

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

    def getRelativePos(self,strand,_start,end,cigar,pos,traceintervalLen):

        #tp = (strand,start,end,cigar,pos,traceintervalLen)

        rel0 = pos - _start - SEQ_TAKE_MARGIN +1
        rel = self.correctCigar(rel0,cigar)
        start = rel
        if start < 2:
            # do not use lower end
            return None

        rel = pos - _start + 6 + SEQ_TAKE_MARGIN +1
        end = self.correctCigar(rel, cigar)
        if end >= traceintervalLen:
            #end = traceintervalLen-1
            return None
        if self.IndelStrict:
            expectedlen = 6 + (SEQ_TAKE_MARGIN*2)
            if abs((end-start) - expectedlen) > 0: # Do not use cross Indel
                #print((end-start))
                return None

        # print(tp)
        #print("start end",cigar,pos,_start,rel0,start,end)
        return start,end


    def calcStartEnd(self,strand,start,end,cigar,pos,offset,traceintervals):

        rp = self.getRelativePos(strand,start,end,cigar,pos,len(traceintervals))
        if rp is None:
            return None
        #print("offset", offset)
        relativeStart, relativeEnd = rp
        # print("relativeStart, relativeEnd",relativeStart, relativeEnd)
        # print(traceintervals)
        if relativeEnd >= len(traceintervals) or relativeStart >= len(traceintervals):
            return None
        #tracestart = traceintervals[relativeStart]
        #traceend = traceintervals[relativeEnd]
        #print(traceintervals)
        #print("sig_start,sig_end",traceintervals[relativeStart],traceintervals[relativeEnd])
        return traceintervals[relativeStart:relativeEnd]


    def getOneRow(self,row,strand, pos,rseq):

        chr = row['chr']
        start = row['start']
        end = row['end']
        cigar = row['cigar']
        offset = row['offset']
        traceintervals = row['traceintervals']
        trace = row['trace']
        #
        #print(traceintervals)
        sted = self.calcStartEnd(strand,start,end,cigar,pos,offset,traceintervals)
        if sted is None:
            return None
        traceItv = sted
        #print(traceItv,sig_start, sig_end)
        if len(traceItv) < 2:
            return None

        tracestart = traceItv[0]
        traceend = traceItv[-1]
        trace = trace.reshape(-1,4)
        trace = trace.T
        trace = trace[::,tracestart:traceend]

        #print(signal)
        #
        binnedTrace = binned(trace,traceItv,DATA_LENGTH_UNIT,rseq)
        return binnedTrace

    def getFormattedData(self, strand, pos,rseq,_takecnt):

        traces = []
        takecnt = 0
        for index, row in self.bufferData.iterrows():

            trace = self.getOneRow(row, strand, pos, rseq)
            if trace is not None:

                traces.append(trace)
                takecnt = takecnt + 1
                #print("takecnt",takecnt,_takecnt)
                if _takecnt > 0 and takecnt == _takecnt:
                    break

        # data = np.array(data)
        # data = data / 256
        return traces,takecnt
