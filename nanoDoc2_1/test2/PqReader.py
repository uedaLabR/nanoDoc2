import pyarrow.parquet as pq
import mappy as mp
import pysam


def getChrom(f):

    return f.split("/")[-3].replace("chrom=","")

def getStrand(f):

    return f.split("/")[-2].replace("strand=","") == "0"

def getFiles(self,path,chrom,strand,findexs):

    sortedfile = []
    for root, subFolder, files in os.walk(path):
        for item in files:
            if item.endswith(".parquet"):
                fileNamePath = str(os.path.join(root, item))
                chromFromPath = getChrom(fileNamePath)
                strandP = getStrand(fileNamePath)

                if (chromFromPath == chrom) and (strand == strandP):
                    sortedfile.append(fileNamePath)

    sortedfile = sorted(sortedfile)
    ret = []
    n = 0
    for f in sortedfile:
        if n in findexs:
            ret.append(f)
        n+=1

    return ret

import os
import pandas as pd
import numpy as np
from scipy import ndimage as ndi
import statistics


def intervalToAbsolute(intervals):
    ret = []
    cnt = 0
    sum = 0
    for n in intervals:

        if cnt == 0:
            ret.append(0)
            sum = sum + n
            cnt += 1
            continue
        # else
        ret.append(sum)
        sum = sum + n

    ret.append(sum)
    return np.array(ret)

class TraceSeqReader:

    def getChrom(self,f):

        return f.split("/")[-2].replace("chrom=", "")

    def getStrand(self,f):

        return f.split("/")[-1].replace("strand=", "")

    def getFiles(self,path,chrom,strand,findexs):

        fidx = 0
        filelist = []
        for root, subFolder, files in os.walk(path):
            for item in files:
                if item.endswith(".parquet"):
                    chrom_f = self.getChrom(root)
                    strand_f = self.getStrand(root)
                    strandok = False
                    if strand and strand_f=="1":
                        strandok = True
                    elif not strand and not (strand_f=="1"):
                        strandok = True

                    if strandok and chrom == chrom_f:
                        filelist.append(root + "/" + item)

        filelist = sorted(filelist)
        filelist2 = []
        idx = 0
        for file in filelist:

            if idx in findexs:
                filelist2.append(file)
            idx+=1

        return filelist2



    def __init__(self,chrom,strand,start, end, path,findexs,maxreads,minreadlen):

        files = self.getFiles(path,chrom,strand,findexs)
        print(files)
        data = None
        for file in files:

            dataadd = pq.read_table(file, columns=['r_st', 'r_en','q_st','q_en','cigar','traceintervals','signal']).to_pandas()
            dataadd['traceintervals'] = dataadd['traceintervals'].apply(intervalToAbsolute)
            if data is None:
                data = dataadd
            else:
                data = pd.concat([data, dataadd])

        self.data = data
        self.strand = strand
        self.bufData = None
        self.minreadlen = minreadlen
        self.start = start
        self.end = end



    def load(self,pos):

        query = 'r_st < ' + str(pos) + ' & r_en > ' + str(pos) + \
                ' & (r_en-r_st) >= ' + str(self.minreadlen)
        self.bufData =  self.data.query(query)
        #

    import pysam
    def correctCigar(self,targetPos,cigar):

        a = pysam.AlignedSegment()
        a.cigarstring = cigar
        refpos = 0
        relpos = 0
        cidx = 0
        for cigaroprator, cigarlen in a.cigar:

            if cigaroprator == 3 and cidx > 0: #N

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

        cidx+=1

        return 0


    def getRelativePos(self, start, end, cigar, pos):

        if self.strand == True:
            return self.getRelativePosP(start, end, cigar, pos)
        else:
            return self.getRelativePosN(start, end, cigar, pos)

    def getRelativePosP(self, _start, end, cigar, pos):

        # tp = (strand,start,end,cigar,pos,traceintervalLen)
        rel0 = pos - _start
        rel = self.correctCigar(rel0, cigar)
        start = rel
        startmargin = 8
        if start < startmargin:
            # do not use lower end
            return None

        rel = pos - _start + 6
        end = self.correctCigar(rel, cigar)
        return start, end

    def getRelativePosN(self, start, end, cigar, pos):

        # tp = (strand,start,end,cigar,pos,traceintervalLen)
        margin = 1
        rel0 = end - pos - margin
        rel = self.correctCigar(rel0, cigar)
        start = rel
        startmargin = 8
        if start < startmargin:
            # do not use lower end
            return None

        rel = end - pos + 6 + margin
        end = self.correctCigar(rel, cigar)
        return start, end

    def calcStartEnd(self,start,end,cigar,pos):

        rp = self.getRelativePos(start,end,cigar,pos)
        if rp is None:
            return None
        # print(rp)
        #print("offset", offset)
        relativeStart, relativeEnd = rp
        if abs(relativeStart-relativeEnd) != 6:
            return None
        return relativeStart,relativeEnd


    def binSignal(self,trimsignal, trimlength, mode=1):

        minlen = 50
        maxlen = trimlength * 1.2
        if len(trimsignal) == trimlength:
            return trimsignal  # not very likely to happen

        elif len(trimsignal) > trimlength:
            # trim from first
            # return downsample(trimsignal, trimlength)
            return None
        # elif len(trimsignal) > trimlength:
        #
        #     return downsample(trimsignal, trimlength)

        elif len(trimsignal) < minlen:

            return None

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


    def getOneRow(self,row, pos):

        start = row['r_st']
        end = row['r_en']
        q_start = row['q_st']
        q_end = row['q_en']

        cigar = row['cigar']
        traceboundary = row['traceintervals']
        signal = row['signal']
        UNIT = 10
        UNITLENGTH = 1024

        sted = self.calcStartEnd(start,end,cigar,pos)
        # print(sted,start,end,cigar,pos)
        if sted is None:
            return None
        relativeStart,relativeEnd = sted
        signal_start = traceboundary[relativeStart]*UNIT
        signal_end = traceboundary[relativeEnd]*UNIT

        subsignal = signal[signal_start:signal_end]

        if len(subsignal) == 0:
            return None

        info = (start, end, q_start, q_end, cigar)

        binsignal = self.binSignal(subsignal,UNITLENGTH)
        if binsignal is None:
            return None

        return  binsignal, info

    def getRowData(self, pos, takecnt=-1):

        reloadUnit = 100
        if self.bufData is None or abs(pos-self.start) % reloadUnit:
            self.load(pos)

        initfilterdata = []
        infos = []
        for index, row in self.bufData.iterrows():

            ret = self.getOneRow(row, pos)
            if ret is not None:
                # print(ret.shape)
                v,i = ret
                initfilterdata.append(v)
                infos.append(i)
            if len(initfilterdata) == takecnt:
                break

        return initfilterdata,infos



            # meta = {'sortkey':'u8', 'read_id':'object','chrom':'object','strand':'u1','readstart':'u4','readseq':'object',
    #         'r_st':'u4','r_en':'u4','q_st':'u4','q_en':'u4','cigar':'object',
    #         'fastq':'object','offset':'u1','traceintervals':'object',
    #         'leftfirst':'u4','leftlast':'u4','traceseq':'object',
    #         'trace':'object','signal':'object','signalboundary':'object',
    #        }