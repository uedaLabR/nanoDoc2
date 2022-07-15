import pyarrow.parquet as pq
import matplotlib.pyplot as plt
import numpy as np
# from tensorflow.keras.utils.training_utils import multi_gpu_model
import tensorflow
import pandas as pd
from tensorflow.keras.callbacks import ModelCheckpoint
import pandas as pd
from Bio import SeqIO
from operator import itemgetter, attrgetter
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq
import numpy as np
import glob
import pysam
import mappy as mp
from numba import jit,u1,i8,f8


def cigarToMap(cigar,strand):

    startadjust = 0
    cigseqlen = 0
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    cigarlist = []
    isFirst = True
    cigarL = a.cigar
    if not strand:
        cigarL.reverse()

    for cigaroprator, cigarlen in a.cigar:

        if isFirst:
            isFirst = False
            if cigaroprator == 2 or cigaroprator == 4:  # Del OR SOFTCLIP
                startadjust = cigarlen

        if  cigaroprator == 0 or cigaroprator == 1 or cigaroprator == 4:
            cigseqlen += cigarlen

        cigarlist.append((cigaroprator, cigarlen))

    return cigarlist,startadjust,cigseqlen


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

def getTraceSeq(trace,signalintervals):

    trace = decode(trace)
    trace = trace.T
    #print(trace)
    sl = []
    seq = ["A","C","G","T"]

    for m in range(1,len(signalintervals)):

        b4 = int(signalintervals[m-1]) // 10
        idx = int(signalintervals[m]) // 10

        if b4 == idx:
            sl.append("N")
            # print("N",m)
            continue

        partialtrace = trace[:,b4:idx]
        su = np.sum(partialtrace, axis=1)
        #print("su",su)
        maxidx = su.argmax()
        s = seq[maxidx]
        sl.append(s)

    return "".join(sl)

def getTraceSeq2(trace,signalintervals):

    trace = decode(trace)
    trace = trace.T
    #print(trace)
    sl = []
    seq = ["A","C","G","T"]

    for m in range(1,len(signalintervals)):

        b4 = int(signalintervals[m-1]) // 10
        idx = int(signalintervals[m]) // 10


        partialtrace = trace[:,b4:idx]
        su = np.sum(partialtrace, axis=1)
        #print("su",su)
        maxidx = su.argmax()
        s = seq[maxidx]
        sl.append(s)

    return "".join(sl)


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

def toBamRecord(reflist,row):

    a = pysam.AlignedSegment()
    a.query_name = row.read_id

    a.flag = 0
    seq = getTraceSeq(row.trace, row.signalboundary)
    # cigar = row.cigar
    cigar = str(len(seq))+"M"
    print(cigar)
    adj = 0
    if row.strand == False:
         a.flag = 16
         seq = mp.revcomp(seq)
         adj = 1


    a.query_sequence = seq

    a.reference_id = reflist.index(row.chr)
    a.cigar,startadjust,cigarseqlen = cigarToMap(cigar,row.strand)
    a.reference_start = row.start + startadjust +adj
    a.mapping_quality = 20

    if a.query_sequence is None or len(a.query_sequence) != cigarseqlen:
        print(row.cigar)
        print(a.query_sequence)
        return None

    return a

def toBamRecord2(reflist,row):

    a = pysam.AlignedSegment()
    a.query_name = row.read_id  +"2"

    a.flag = 0
    seq = getTraceSeq2(row.trace, row.signalboundary)
    cigar = row.cigar

    adj = 0
    if row.strand == False:
         a.flag = 16
         seq = mp.revcomp(seq)
         adj = 1


    a.query_sequence = seq

    a.reference_id = reflist.index(row.chr)
    a.cigar,startadjust,cigarseqlen = cigarToMap(cigar,row.strand)
    a.reference_start = row.start + startadjust +adj
    a.mapping_quality = 20

    if a.query_sequence is None or len(a.query_sequence) != cigarseqlen:
        print(row.cigar)
        print(a.query_sequence)
        return None

    return a


def toBam(ref,path,pathout):

    sortedfile = sorted(glob.glob(path + "/*.pq"))
    # print(sortedfile)
    pluslist = []
    minuslist = []
    for f in sortedfile:
       if "_-1_" in f:
           minuslist.append(f)
       else:
           pluslist.append(f)
           break

    pluslist.sort(key=lambda x: int(x.split('_')[-1].replace(".pq","")))
    minuslist.sort(key=lambda x: int(x.split('_')[-1].replace(".pq","")))
    flist = pluslist
    flist.extend(minuslist)

    sqlist = []
    reflist = []
    for record in SeqIO.parse(ref, 'fasta'):
        sqlist.append({'LN': len(record), 'SN': record.id})
        reflist.append(record.id)

    header = {'HD': {'VN': '1.0'},
              'SQ': sqlist}

    out_bamfile = pysam.AlignmentFile(pathout, "wb", header=header)

    for f in flist:

        print(f)
        parquet_file = pq.ParquetFile(f)
        columns = ['read_no', 'read_id','chr', 'strand', 'start', 'end', 'cigar','genome','offset', 'trace','signal','signalboundary']
        df = pq.read_table(f, columns=columns).to_pandas()
        df['signalboundary'] = df['signalboundary'].apply(intervalToAbsolute)

        for index, row in df.iterrows():

            bamrecord = toBamRecord(reflist,row)
            if bamrecord is not None:
                out_bamfile.write(bamrecord)
            bamrecord2 = toBamRecord2(reflist,row)
            if bamrecord2 is not None:
                out_bamfile.write(bamrecord2)

    out_bamfile.close()



from Bio import SeqIO
if __name__ == "__main__":


    ref = "/data/nanopore/reference/Curlcake.fa"
    path = '/data/nanoDoc2_2/varidate/CurlcakeIVT'
    pathout = "/data/nanoDoc2_2/Curlcake.bam"
    toBam(ref,path,pathout)
