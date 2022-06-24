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


@jit
def decode16bit(a_trace):

    a = (a_trace & 0b1111000000000000) >> 12
    c = (a_trace & 0b0000111100000000) >> 8
    g = (a_trace & 0b0000000011110000) >> 4
    t = (a_trace & 0b0000000000001111)
    #
    return (a,c,g,t)

import numpy as np
def decode(trace):

    tp = list(map(decode16bit, trace))
    ar = np.array(tp)
    return ar


from scipy import ndimage as ndi
import statistics
def eachScore(nuc,atrace):


    if nuc == 'A':
        nucidx = 0
    elif nuc == 'C':
        nucidx = 1
    elif nuc == 'G':
        nucidx = 2
    else:
        nucidx = 3

    #print("atrace",atrace)
    atrace = decode(atrace)
    #print("decode",atrace)
    if len(atrace) > 1:
        rowsum =np.sum(atrace, axis=0)
    else:
        rowsum = atrace[0]

    #print("rowsum",rowsum)
    sum = np.sum(rowsum)
    matchsum = rowsum[nucidx]
    return matchsum/sum

import nanoDoc2_1.preprocess.WriteToFile as wr
def getScore(subtraces,seq):

    score = 0
    upto = min(len(seq),len(subtraces))
    #print(subtraces)
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
        idx = n-1
        if score > maxscore:
            maxidx = idx

    #print("maxidx",maxidx)
    return maxidx

import scipy.signal as scisignal
from scipy.interpolate import interp1d
def downsample(array, npts):

    interpolated = interp1d(np.arange(len(array)), array, axis = 0, fill_value = 'extrapolate')
    downsampled = interpolated(np.linspace(0, len(array), npts))
    #downsampled = scisignal.resample(array, npts)
    return downsampled


def binSignal(trimsignal, trimlength, mode=1):

    if len(trimsignal) == trimlength:
        return trimsignal  # not very likely to happen

    if len(trimsignal) > trimlength:
        # trim from first
        return downsample(trimsignal, trimlength)

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



        # rand1 = np.random.rand(lefthalf) * sigma
        # rand1 = rand1 + med
        leftlen = trimlength - siglen - lefthalf
        # rand2 = np.random.rand(leftlen) * sigma
        # rand2 = rand2 + med
        # #
        # ret = np.concatenate([rand1, trimsignal, rand2])

        ret = np.concatenate([np.zeros(lefthalf), trimsignal, np.zeros(leftlen)])

        return ret

def binTrace(trace, trimlength):

    x = len(trace)
    if x == trimlength:
        return trace
    elif x > trimlength:
        half = (x-trimlength)//2-1
        return trace[half:half+trimlength]
    else:
        half = (trimlength-x)//2 - 1
        left = trimlength - half- x
        if half < 0:
            half = 0
        if left < 0:
            left  = 0
        return np.pad(trace, (half,left), 'constant')


DATA_LENGTH = 768
DATA_LENGTH_Trace = 80
def binned(trimtrace, traceItv,trimUnitLength,rseq,signal):

    bintrace = binTrace(trimtrace,DATA_LENGTH_Trace)
    if len(signal) == 0:
        return None
    binsignal = binSignal(signal, DATA_LENGTH)

    return bintrace,binsignal


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



    def __init__(self, path,ref,minreadlen,strand,start=0, end=0, maxreads = 500,IndelStrict = False):

        self.IndelStrict = IndelStrict
        self.path = path
        self.batch = None
        self.strand = strand
        print("ref",ref)
        self.a = mp.Aligner(ref)
        self.maxreads = int(maxreads * 1.2)# take little more sample since some reads disqualify
        self.maxreads_org = maxreads
        self.minreadlen = minreadlen
        self.bufferData = None
        self.loadchr = None
        self.firstbinsizeList = [15,30,100,200]
        self.takemarginb4 = 2
        self.takemargin = 2
        self.lastloadpos = 0
        self.firstloadpos = 0
        self.binsize = 5000
        self.samplelen = 4000
        self.start = start
        self.end = end


        #make index from parquet meta data
        indexlist = []
        fileidx = 0
        sortedfile = self.getFilePathList(path)

        print("sortedfile",sortedfile)
        for file in sortedfile:

            parquet_file = pq.ParquetFile(file)

            chrInfo = parquet_file.metadata.row_group(0).column(2).statistics.min
            strandInfo = parquet_file.metadata.row_group(0).column(3).statistics.min
            startInfo = parquet_file.metadata.row_group(0).column(4).statistics.min
            endInfo = parquet_file.metadata.row_group(0).column(5).statistics.max
            indexlist.append((fileidx,chrInfo,strandInfo,startInfo,endInfo))
            fileidx +=1

        self.indexdf = pd.DataFrame(indexlist,
                          columns=['fileidx','chr','strand','start','end'])
        #print(self.indexdf)

    def getDepth(self,chr, pos, strand):


        query = 'start <= ' + str(pos-self.takemarginb4) + ' & end >= ' + str(pos+self.takemargin) + \
                ' & chr == "' + chr + '" & strand == ' + str(strand) + '' + \
                ' & (end-start) >= ' + str(self.minreadlen)

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
                #print(filepath)

        #print(indexdata)
        depth = indexdata.query('start <=' + str(pos - 10) + ' & end >=' + str(pos + 10))['start'].count()
        #depth = indexdata.query('start <=' + str(pos-10) + ' & end >=' + str(pos+10))['start'].count()
        return depth


    def randomsample(self, pos, data, ntake, indexes):

        if indexes is None:

            #print('init read')
            if data is None:
                return None
            datainpos = data.query('start <=' + str(pos-self.takemarginb4) + ' & end >=' + str(pos+self.takemargin) +" & (end-start) > "+str(self.minreadlen))

            if datainpos is None or len(datainpos) ==0:
                return None
            leno = len(datainpos)
            if ntake < len(datainpos):
                datainpos = datainpos.sample(n=ntake)
                #datainpos = datainpos.head(ntake)
                print("random sampling 1",pos,leno,ntake)

            return datainpos

        else:

            dataposprev = indexes.query('start <=' + str(pos-self.takemarginb4) + ' & end >=' + str(pos+self.takemargin)+" & (end-start) > "+str(self.minreadlen))
            if len(dataposprev) <= ntake:
                return None

            else:
                datainpos = data.query('start <=' + str(pos-self.takemarginb4) + ' & end >=' + str(pos+self.takemargin)+" & (end-start) > "+str(self.minreadlen))

                df_alreadyhave = dataposprev['read_no']
                datainpos = datainpos[~datainpos.read_no.isin(df_alreadyhave)]
                cnt = ntake - len(df_alreadyhave)
                if (cnt > 0) and (cnt <= len(datainpos)):
                    datainpos = datainpos.sample(n=cnt)
                    #datainpos = datainpos.head(cnt)
                    return datainpos

                return None



    def load(self,chr, pos, strand):

        print("loding data")
        self.loadchr = chr

        query = 'start <= ' + str(pos-self.takemarginb4) + ' & end >= ' + str(pos+self.takemargin) + \
                ' & chr == "' + chr + '" & strand == ' + str(strand) + '' + \
                ' & (end-start) >= ' + str(self.minreadlen)

        print(query)
        print(self.indexdf)
        pqfiles = self.indexdf.query(query)
        #
        sortedfile = self.getFilePathList(self.path)
        #print(sortedfile)

        #
        indexdata = None
        for index, row in pqfiles.iterrows():

            fileidx = row['fileidx']
            filepath = sortedfile[fileidx]
            if indexdata is None:

                indexdata = pq.read_table(filepath, columns=['read_no','read_id', 'chr', 'strand', 'start', 'end']).to_pandas()
                indexdata['fileidx'] = fileidx

            else:
                dataadd = pq.read_table(filepath, columns=['read_no','read_id', 'chr', 'strand', 'start', 'end']).to_pandas()
                dataadd['fileidx'] = fileidx
                indexdata = pd.concat([indexdata, dataadd])

        # get reads ids for bufferted position (start,end)
        start = pos
        end = pos +  self.samplelen

        if not strand:
            start = pos - self.samplelen
            end = pos

        readsIndex = None
        ntake = self.maxreads
        if self.firstloadpos == 0:
            self.firstloadpos = pos
        self.lastloadpos = pos

        rlist = list(range(start,end))
        if not strand:
            rlist.reverse()

        for pos2 in rlist:

            addIndex = self.randomsample(pos2, indexdata, ntake, readsIndex)
            if readsIndex is None:
                readsIndex = addIndex
            elif addIndex is not None:

                print("addIndex")
                print(pos2, addIndex)
                readsIndex = pd.concat([readsIndex, addIndex])


        if readsIndex is None:
            return None

        # read data with row signal
        fileindexes = readsIndex['fileidx'].unique()
        dataWithRow = None
        for findex in fileindexes:

            filepath = sortedfile[findex]
            readids = readsIndex['read_no']
            filterlist = []
            filterTp = ('read_no', 'in', readids)
            filterlist.append(filterTp)
            columns = ['read_no', 'read_id','chr', 'strand', 'start', 'end', 'cigar','genome','offset', 'traceintervals','trace','signal']
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

        rellastload = pos-self.lastloadpos
        relfirstload = pos-self.firstloadpos
        if not strand:
            rellastload = self.lastloadpos - pos
            relfirstload = self.firstloadpos - pos

        firstload = (relfirstload in self.firstbinsizeList)
        nonlast = (self.end-pos) > 1000
        loadpos = (rellastload %self.binsize == 0) and rellastload > 0 and nonlast
        #print(relfirstload,rellastload,firstload,loadpos)
        if ((self.bufferData is None) or (chr != self.loadchr) or firstload or loadpos):
            self.bufferData = None
            self.load(chr, pos, strand)

        traces,signals,sampledlen,infos = self.getFormattedData(strand, pos,rseq,takecnt)
        if sampledlen > self.maxreads_org:
            sampledlen = self.maxreads_org
            signals = signals[0:self.maxreads_org]

        if sampledlen < takecnt and takecnt > 0:
            depth = self.getDepth(chr, pos, strand)
            if depth > sampledlen + 50:
                #load and sample again since not enough sampling for this region
                print("load data",depth,sampledlen,takecnt)
                self.load(chr, pos, strand)
                traces,signals,sampledlen,infos = self.getFormattedData(strand, pos,rseq, takecnt)

        if takecnt == -1:
            depth = self.getDepth(chr, pos, strand)
            if sampledlen < depth and nonlast:
                self.load(chr, pos, strand)
                traces,signals,sampledlen,infos = self.getFormattedData(strand, pos,rseq, takecnt)

        return traces,signals,sampledlen,infos


    def getRowSequence(self, chr, strand, pos,takecnt=-1):

        binsizehalf = self.binsize //2
        if ((self.bufferData is None) or (chr != self.loadchr) or ((pos % binsizehalf) == 0)):
            self.load(chr, pos, strand)

        return self.bufferData.iterrows()

    import pysam
    def correctCigar(self,targetPos,cigar):

        a = pysam.AlignedSegment()
        a.cigarstring = cigar
        refpos = 0
        relpos = 0
        for cigaroprator, cigarlen in a.cigar:

            if cigaroprator == 3: #N

                refpos = refpos + cigarlen

            elif cigaroprator == 0 or cigaroprator == 4:  # match or S softclip was not correted so treat as M

                if refpos + cigarlen > targetPos:
                    return relpos + (targetPos - refpos)

                relpos = relpos + cigarlen
                refpos = refpos + cigarlen

            elif cigaroprator == 2:  # Del

                refpos = refpos + cigarlen

            elif cigaroprator == 1 :  # Ins

                if relpos == 0:
                    if targetPos <= cigarlen:
                        return 0

                relpos = relpos + cigarlen

        return 0

    def getRelativePos(self, strand, _start, end, cigar, pos, traceintervalLen):

        if strand:
            return self.getRelativePosP(strand,_start,end,cigar,pos,traceintervalLen)
        else:
            return self.getRelativePosN(strand, _start, end, cigar, pos, traceintervalLen)

    def getRelativePosP(self,strand,_start,end,cigar,pos,traceintervalLen):

        #tp = (strand,start,end,cigar,pos,traceintervalLen)
        margin = 1
        rel0 = pos - _start - margin
        rel = self.correctCigar(rel0,cigar)
        start = rel
        startmargin = 8
        endmargin = 10
        if start < startmargin:
            # do not use lower end
            return None

        rel = pos - _start + 7 +  margin
        end = self.correctCigar(rel, cigar)
        #do not take last 5
        if end == 0 or end >= traceintervalLen-endmargin:
            #end = traceintervalLen-1
            return None
        if self.IndelStrict:
            expectedlen = 9
            if abs((end-start) - expectedlen) > 0: # Do not use cross Indel
                return None

        # print(tp)
        #print("start end",cigar,pos,_start,rel0,start,end)
        return start,end

    def getRelativePosN(self,strand,_start,end,cigar,pos,traceintervalLen):

        #tp = (strand,start,end,cigar,pos,traceintervalLen)
        margin = 1
        rel0 = end - pos - margin
        rel = self.correctCigar(rel0,cigar)
        start = rel
        startmargin = 8
        endmargin = 10
        if start < startmargin:
            # do not use lower end
            return None

        rel = end - pos + 7 +  margin
        end = self.correctCigar(rel, cigar)
        #do not take last 5
        if end == 0 or end >= traceintervalLen-endmargin:
            #end = traceintervalLen-1
            return None
        if self.IndelStrict:
            expectedlen = 9
            if abs((end-start) - expectedlen) > 0: # Do not use cross Indel
                return None

        # print(tp)
        #print("start end",cigar,pos,_start,rel0,start,end)
        return start,end

    def getScore(self,subtraces, seq):

        score = 0
        upto = min(len(seq), len(subtraces))
        for n in range(upto):
            nuc = seq[n]
            atrace = subtraces[n]
            ascore = eachScore(nuc, atrace)
            score += ascore
        return score


    def analyzeIdxShift(self,subtraces, trace,subtraceIntv,rseq):


        maxidx = 0
        maxscore = 0
        for n in range(6):
            seq = rseq[n:n + 10]
            score = self.getScore(subtraces, seq)
            idx = n - 3
            if score > maxscore:
                maxidx = idx

        return maxidx

    def calcStartEnd(self,strand,start,end,cigar,pos,offset,traceintervals,trace,rseq):

        rp = self.getRelativePos(strand,start,end,cigar,pos,len(traceintervals))
        if rp is None:
            return None
        #print("offset", offset)
        relativeStart, relativeEnd = rp
        if relativeEnd >= len(traceintervals) or relativeStart >= len(traceintervals):
            return None
        subtraceIntv = traceintervals[relativeStart:relativeEnd]
        return subtraceIntv
        # shift = self.ananlyzeshift(trace,subtraceIntv,rseq)
        # return traceintervals[relativeStart+shift:relativeEnd+shift]


    def getOneRow(self,row,strand, pos,rseq):

        id = row['read_id']
        chr = row['chr']
        start = row['start']
        end = row['end']
        cigar = row['cigar']
        offset = row['offset']
        traceintervals = row['traceintervals']
        trace = row['trace']
        signal = row['signal']
        #
        #print(traceintervals)
        sted = self.calcStartEnd(strand,start,end,cigar,pos,offset,traceintervals,trace,rseq)
        if sted is None:
            return None
        traceItv = sted
        #print(traceItv,sig_start, sig_end)
        if len(traceItv) < 2:
            return None
        SignalUNIT = 10
        #print(traceItv)
        tracestart = traceItv[0]
        traceend = traceItv[-2]
        trace_t = trace[tracestart:traceend]
        signal_t = signal[tracestart*SignalUNIT:traceend*SignalUNIT]


        #
        ret = binned(trace_t,traceItv,DATA_LENGTH_UNIT,rseq,signal_t)
        if ret is None:
            return None
        info = chr+":"+str(start)+"-"+str(end)+" "+cigar
        bintrace,binsignal = ret
        return bintrace,binsignal,info

    def getFormattedData(self, strand, pos,rseq,_takecnt):

        untrimsignals = []
        untrimtraces = []
        signals = []
        traces = []
        infos = []
        ids = []
        takecnt = 0
        for index, row in self.bufferData.iterrows():

            ret = self.getOneRow(row, strand, pos, rseq)

            if ret is not None:

                trace,signal,info  = ret
                traces.append(trace)
                signals.append(signal)
                infos.append(info)

                takecnt = takecnt + 1
                #print("takecnt",takecnt,_takecnt)
                if _takecnt > 0 and takecnt == _takecnt:
                    break

        return  traces,signals,takecnt,infos

