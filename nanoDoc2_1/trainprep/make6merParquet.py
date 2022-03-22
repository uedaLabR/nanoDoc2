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

from nanoDoc2_1.utils.PqFileReader import PqReader
import sys
import random

def makeSamplePlan(refs,pqs,output_file,takeCnt):

    fivemerDict={}
    recordL=[]
    cnt = 0
    for ref in refs:
        records = SeqIO.parse(ref, 'fasta')
        fr = PqReader(pqs[cnt], ref,3000,IndelStrict=True)

        for record in records:
           print(record.name)
           print(record.seq)
           seq = record.seq

           for n in range(30, len(seq)-30):
           #for n in range(30, 100):

              smer = seq[n:n+6]
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
    cntloop = 0
    for p in posList:

        fileidx,chr,pos,depth,takecnt,fmer = p

        #
        if pfidx != fileidx:
            path = pqs[fileidx]
            ref = refs[fileidx]
            print("init reader")
            fr = PqReader(path,ref, 3000)

        traces,signals,sampledlen = fr.getRowData(chr, True, pos,takecnt=takecnt)
        print(p,len(signals),takecnt)
        pfidx = fileidx
        if len(signals) > 0:
            if fmer in datadict:
                datadict[fmer].extend(signals)
            else:
                d = []
                d.extend(signals)
                datadict[fmer] = d

    #write to parquet
    keys = datadict.keys()
    keys = sorted(keys)
    dataf = []
    keyidx = 0
    for key in keys:

        d = datadict[key]
        random.shuffle(d)
        for da in d:
            da = da.flatten()
            tp = (keyidx,key,da)
            dataf.append(tp)
            keyidx += 1


    pschema = schema(
        [
            ('flg', uint32()),
            ('fmer', string()),
            ('signal', list_(float32()))
        ]
    )
    df = pd.DataFrame(dataf,
                      columns=['flg','fmer', 'signal'])
    pd.set_option('display.max_rows', None)


    pyarrow_table = Table.from_pandas(df, pschema)
    pq.write_table(
        pyarrow_table,
        output_file,
        row_group_size=4000,
        compression='snappy',
        flavor=['spark'],
    )

        # sidx = 0
        # sb4 = None
        #
        # f = open(path_w, mode='w')
        # for d in posList:
        #     f.write(",".join(map(str, d)) + "\n")
        # f.close()





    # outdir = ""
    # args = sys.argv
    # path = args[1]
    # print(path)

