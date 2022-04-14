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


def cigarToMap(cigar):

    startadjust = 0
    cigseqlen = 0
    a = pysam.AlignedSegment()
    a.cigarstring = cigar
    cigarlist = []
    isFirst = True
    for cigaroprator, cigarlen in a.cigar:

        if isFirst:
            isFirst = False
            if cigaroprator == 2 or cigaroprator == 4:  # Del OR SOFTCLIP
                startadjust = cigarlen

        if  cigaroprator == 0 or cigaroprator == 1 or cigaroprator == 4:
            cigseqlen += cigarlen

        cigarlist.append((cigaroprator, cigarlen))

    return cigarlist,startadjust,cigseqlen


def getTraceSeq(traceintervals,trace):

    #print("ti",traceintervals)
    trace = trace.reshape(-1, 4)
    trace = trace.T
    #print("trace", trace)
    sl = []
    seq = ["A","C","G","T"]
    for m in range(1,len(traceintervals)):

        b4 = traceintervals[m-1]
        idx = traceintervals[m]
        partialtrace = trace[:,b4:idx]
        #print("partialtrace",partialtrace)
        su = np.sum(partialtrace, axis=1)
        #print(su)
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
    # columns = ['read_id', 'chr', 'strand', 'start', 'end', 'cigar', 'genome', 'offset', 'traceintervals', 'trace']
    a.query_name = row.read_id
    a.query_sequence = getTraceSeq(row.traceintervals,row.trace)
    a.flag = 99
    a.reference_id = reflist.index(row.chr)
    a.cigar,startadjust,cigarseqlen = cigarToMap(row.cigar)
    a.reference_start = row.start + startadjust
    a.mapping_quality = 20
    # print(row.cigar)
    # print(a.query_sequence)
    # print(len(a.query_sequence),cigarseqlen)
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

    #print(pluslist)
    #print(minuslist)
    for f in flist:

        print(f)
        parquet_file = pq.ParquetFile(f)
        columns = ['read_id', 'chr', 'strand', 'start', 'end', 'cigar', 'genome', 'offset', 'traceintervals', 'trace']
        df = pq.read_table(f, columns=columns).to_pandas()
        df['traceintervals'] = df['traceintervals'].apply(intervalToAbsolute)
        for index, row in df.iterrows():

            bamrecord = toBamRecord(reflist,row)
            if bamrecord is not None:
                out_bamfile.write(bamrecord)

    out_bamfile.close()



from Bio import SeqIO
if __name__ == "__main__":

    ref = "/data/nanopore/reference/NC000913.fa"
    path = '/data/nanopore/nanoDoc2/1623_ivt'
    pathout = "/data/nanopore/nanoDoc2/1623_ivt.bam"
    toBam(ref,path,pathout)
