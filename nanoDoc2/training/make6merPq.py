import pyarrow.parquet as pq
import pandas as pd
from Bio import SeqIO

class Counter:

    cnt = 1
    s = set()
    l = []

    def __init__(self, v,name,n,depth,idxs):
       self.s = set()
       self.l = []
       self.s.add(v)
       self.l.append((name,n,depth,idxs))

    def inc(self,v,name,n,depth,idxs):
       self.cnt = self.cnt + 1
       self.s.add(v)
       self.l.append((name, n, depth,idxs))

    def getS(self):
        return self.s
    def getCnt(self):
        return self.cnt
    def getList(self):
        return self.l

from nanoDoc2.utils.PqFileReader import PqReader
import sys
def makeSamplePlan(refs,pqs,path_w):

    sixmerDict={}
    recordL=[]
    cnt = 0
    for ref in refs:
        records = SeqIO.parse(ref, 'fasta')
        fr = PqReader(pqs[cnt], 100)

        for record in records:
           print(record.name)
           print(record.seq)
           seq = record.seq

           for n in range(10, len(seq)-16):

              smer = seq[n:n+6]
              b4 = seq[n-1]
              after = seq[n+7]
              v = b4+after
              depth = fr.getDepth(record.name,n,True)
              print(n,smer,depth)
              #
              # if sixmerDict.get(smer) == None:
              #     sixmerDict[smer] = Counter(v,record.name,n,depth)
              # else:
              #     sixmerDict[smer].inc(v,record.name,n,depth)

        cnt +=1

if __name__ == "__main__":

    ref1 = "/data/nanopore/reference/Curlcake.fa"
    ref2 = "/data/nanopore/reference/Cov2_Korea.fa"

    path_w = ""
    refs= [ref1]
    pq1 = "/data/nanopore/nanoDoc2/testCurlcakeIVT"
    pqs = [pq1]
    #
    makeSamplePlan(refs,pqs, path_w)

    # outdir = ""
    # args = sys.argv
    # path = args[1]
    # print(path)

