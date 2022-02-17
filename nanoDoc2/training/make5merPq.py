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

from nanoDoc2.utils.PqFileReader import PqReader
import sys
import random

def makeSamplePlan(refs,pqs,output_file,takeCnt):

    fivemerDict={}
    recordL=[]
    cnt = 0
    for ref in refs:
        records = SeqIO.parse(ref, 'fasta')
        fr = PqReader(pqs[cnt], 4000)

        for record in records:
           print(record.name)
           print(record.seq)
           seq = record.seq

           for n in range(30, len(seq)-30):
           #for n in range(30, 100):

              smer = seq[n:n+5]
              b4 = seq[n-1]
              after = seq[n+6]
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
            print("init reader")
            fr = PqReader(path, 4000)

        data, cnt = fr.getRowData(chr, True, pos,takecnt=takecnt)
        print(p,len(data),cnt)
        pfidx = fileidx
        if len(data) > 0:
            if fmer in datadict:
                datadict[fmer].extend(data)
            else:
                d = []
                d.extend(data)
                datadict[fmer] = d

    #write to parquet
    keys = datadict.keys()
    keys = sorted(keys)
    dataf = []
    keyidx = 0
    for key in keys:

        d = datadict[key]
        random.shuffle(d)
        d = np.ravel(d)
        tp = (keyidx,key,d)
        dataf.append(tp)
        keyidx +=1

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

if __name__ == "__main__":

    ref1 = "/data/nanopore/reference/Curlcake.fa"
    ref2 = "/data/nanopore/reference/Cov2_Korea.fa"

    path_w = "/data/nanopore/nanoDoc2/5000each.pq"
    refs= [ref1]
    pq1 = "/data/nanopore/nanoDoc2/testCurlcakeIVT"
    pq2 = "/data/nanopore/nanoDoc2/testSARSCOV2"
    pqs = [pq1]
    #

    # fr = PqReader(pq1, 4000)
    # chr = "cc6m_2244_t7_ecorv"
    # data = fr.getRowData(chr, True, 30, 50)

    takeCnt = 5000
    makeSamplePlan(refs,pqs, path_w,takeCnt)



    # outdir = ""
    # args = sys.argv
    # path = args[1]
    # print(path)

