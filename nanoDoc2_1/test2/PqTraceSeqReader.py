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

            dataadd = pq.read_table(file, columns=['r_st', 'r_en','q_st','q_en','cigar','traceseq']).to_pandas()
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

    def getOneRow(self,row, pos):

        start = row['r_st']
        end = row['r_en']
        q_start = row['q_st']
        q_end = row['q_en']

        cigar = row['cigar']
        traceseq = row['traceseq']
        traceseq = traceseq.reshape(-1, 6)

        sted = self.calcStartEnd(start,end,cigar,pos)
        # print(sted,start,end,cigar,pos)
        if sted is None:
            return None
        relativeStart,relativeEnd = sted
        # print(relativeStart,relativeEnd)
        # print(len(traceseq),traceseq.shape)
        trace_s = traceseq[relativeStart:relativeEnd,:].copy()
        trace_s[:,5] = trace_s[:,5]/256
        info = (start,end,q_start,q_end,cigar)
        return  trace_s, info

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