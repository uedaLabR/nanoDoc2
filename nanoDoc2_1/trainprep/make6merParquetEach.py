import os

import pyarrow.parquet as pq
import pandas as pd
from Bio import SeqIO
from operator import itemgetter, attrgetter
from pyarrow import list_, bool_ ,int8,uint8, uint16, uint32, int64, float64, float32, bool_, date32, decimal128, timestamp, string, Table, schema, parquet as pq
import numpy as np

class Counter:

    cnt = 1
    s = set()
    l = []

    def __init__(self, v,name,n,depth,fidx):
       self.s = set()
       self.l = []
       self.s.add(v)
       self.l.append((name,n,depth,fidx))

    def inc(self,v,name,n,depth,fidx):
       self.cnt = self.cnt + 1
       self.s.add(v)
       self.l.append((name, n, depth,fidx))

    def getS(self):
        return self.s
    def getCnt(self):
        return self.cnt
    def getList(self):
        return self.l

from nanoDoc2_1.utils.PqFile6merReader import PqReader
import sys
import random
import mappy as mp
from numba import jit

def toTuple(tracel, signall,key):
    # print(signall)
    dataf = []
    keyidx = 0
    for n in range(len(tracel)):
        trace = tracel[n]
        signal = signall[n]
        signal = signal.flatten()
        trace = trace.flatten()
        tp = (keyidx, key, trace, signal)
        dataf.append(tp)
        keyidx += 1

    return dataf

def toData(datadict):

    keys = datadict.keys()
    outdict = {}
    keyidx = 0
    for key in keys:

        (tracel, signall) = datadict[key]
        dataf = toTuple(tracel, signall, key)
        outdict[key] = dataf

    return outdict

def bufferout(output_file, datadict, writecnt):

    keys = datadict.keys()

    keyidx = 0
    pschema = schema(
        [
            ('flg', uint32()),
            ('fmer', string()),
            ('trace', list_(uint16())),
            ('signal', list_(float32()))
        ]
    )

    print("buffer writing", len(keys))
    outdict = toData(datadict)
    keys = outdict.keys()
    for key in keys:

        dataf = outdict[key]
        df = pd.DataFrame(dataf,
                      columns=['flg','fmer','trace','signal'])

        pd.set_option('display.max_rows', None)
        pyarrow_table = Table.from_pandas(df, pschema)

        outd = output_file +"/"+key
        if not os.path.exists(outd):
            os.mkdir(outd)
        outf = outd+"/"+str(writecnt)+".pq"

        pq.write_table(
            pyarrow_table,
            outf,
            row_group_size=4000,
            compression='snappy',
            flavor=['spark'],
        )

import glob
def mergePq(path):

    dirs = os.listdir(path=path)

    for dir in dirs:

        fout = path + "/" + dir + ".pq"
        f = path + "/" + dir
        files = glob.glob(f + "/*.pq")
        dftotal = None
        for file in files:

            df = pq.read_table(file).to_pandas()
            if dftotal is None:
                dftotal = df
            else:
                dftotal = pd.concat([dftotal, df])

        pschema = schema(
            [
                ('flg', uint32()),
                ('fmer', string()),
                ('trace', list_(uint16())),
                ('signal', list_(float32()))
            ]
        )

        pyarrow_table = Table.from_pandas(dftotal, pschema)
        pq.write_table(
            pyarrow_table,
            fout,
            row_group_size=4000,
            compression='snappy',
            flavor=['spark'],
        )

def makeSamplePlan(refs,pqs,output_file,takeCnt):

    fivemerDict={}
    recordL=[]
    cnt = 0
    strand = True
    for ref in refs:

        a = mp.Aligner(ref)
        records = SeqIO.parse(ref, 'fasta')
        fr = PqReader(pqs[cnt], ref,400,strand,maxreads = takeCnt,IndelStrict=True)

        for record in records:
           print(record.name)
           print(record.seq)
           seq = record.seq

           for n in range(30, len(seq)-30):
           #for n in range(30, 100):

              smer = seq[n:n+6]
              rseq = a.seq(record.name, start=n, end=n + 6)
              print(smer,rseq)

              b4 = seq[n-1]
              after = seq[n+7]
              v = b4+after
              depth = fr.getDepth(record.name,n,True)
              print(record.name,n,smer,depth)
              #
              if fivemerDict.get(smer) == None:
                  fivemerDict[smer] = Counter(v,record.name,n,depth,cnt)
              else:
                  fivemerDict[smer].inc(v,record.name,n,depth,cnt)

        cnt +=1

    posList = []
    rets = sorted(fivemerDict.items(), key=lambda x: x[0])
    for key,r in rets:
        ls = r.getList()
        llen = len(ls)
        divn = takeCnt // llen
        modn = takeCnt % llen
        #
        dlist = []
        cnt = 0
        for d in ls:
            if cnt == 0:
                dlist.append((d[3],d[0], d[1], d[2], (divn + modn), str(key)))
            else:
                dlist.append((d[3],d[0], d[1], d[2], divn, str(key)))
            cnt = cnt + 1

        posList.extend(dlist)

    posList = sorted(posList, key=itemgetter(0, 1, 2))
    #start reading files
    pfidx = -1

    print("start reading row")
    datadict = {}
    poscnt = 0
    writecnt = 0
    for p in posList:
        poscnt +=1
        fileidx,chr,pos,depth,takecnt,fmer = p

        #
        if pfidx != fileidx:
            path = pqs[fileidx]
            ref = refs[fileidx]
            print("init reader")
            fr = PqReader(path, ref, 400, strand, maxreads=takeCnt,IndelStrict=True)


        traces,signals,sampledlen,infos = fr.getRowData(chr, True, pos,takecnt=takecnt)
        print(p,len(signals),takecnt)

        pfidx = fileidx
        if len(signals) > 0:
            if fmer in datadict:

                (tracel,signall) = datadict[fmer]
                tracel.extend(traces)
                signall.extend(signals)

            else:
                tracel = []
                signall = []
                tracel.extend(traces)
                signall.extend(signals)
                datadict[fmer] = (tracel,signall)

            if poscnt %1000 == 0:
                writecnt+=1
                bufferout(output_file,datadict,writecnt)
                datadict = {}

    writecnt += 1
    bufferout(output_file, datadict, writecnt)
    mergePq(output_file)

    #write to parquet
    # keys = datadict.keys()
    # keys = sorted(keys)
    # print(keys)
    # print(len(keys))
    # dataf = []
    # keyidx = 0
    # for key in keys:
    #
    #     (tracel, signall) = datadict[key]
    #     #print(signall)
    #     print("data len", len(tracel))
    #
    #     for n in range(len(tracel)):
    #
    #         trace = tracel[n]
    #         signal = signall[n]
    #
    #         signal = signal.flatten()
    #         trace = trace.flatten()
    #         tp = (keyidx,key,trace,signal)
    #         dataf.append(tp)
    #         keyidx += 1
    #
    #
    #     pschema = schema(
    #         [
    #             ('flg', uint32()),
    #             ('fmer', string()),
    #             ('trace', list_(uint16())),
    #             ('signal', list_(float32()))
    #         ]
    #     )
    #     df = pd.DataFrame(dataf,
    #                       columns=['flg','fmer','trace','signal'])
    #     pd.set_option('display.max_rows', None)
    #
    #     outf = output_file +"/"+ str(key) +".pw"
    #     pyarrow_table = Table.from_pandas(df, pschema)
    #     pq.write_table(
    #         pyarrow_table,
    #         outf,
    #         row_group_size=4000,
    #         compression='snappy',
    #         flavor=['spark'],
    #     )
    #
    #     # sidx = 0
    #     # sb4 = None
    #     #
    #     # f = open(path_w, mode='w')
    #     # for d in posList:
    #     #     f.write(",".join(map(str, d)) + "\n")
    #     # f.close()





    # outdir = ""
    # args = sys.argv
    # path = args[1]
    # print(path)

